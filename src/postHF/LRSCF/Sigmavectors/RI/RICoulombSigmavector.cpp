/**
 * @file RICoulombSigmavector.cpp
 *
 * @date Dec 07, 2018
 * @author Michael Boeckers, Johannes Toelle
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h"

/* Include Serenity Internal Headers */
#include "basis/BasisFunctionMapper.h" //Construct joined fitting basis.
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RICoulombSigmavector<SCFMode>::RICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                    std::vector<Eigen::MatrixXd> b)
  : Sigmavector<SCFMode>(lrscf, b) {
}

template<Options::SCF_MODES SCFMode>
RICoulombSigmavector<SCFMode>::RICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : Sigmavector<SCFMode>(lrscf) {
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
RICoulombSigmavector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector:   RI-J");

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Thread safety.
  unsigned nThreads = omp_get_max_threads();
  std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> fock_threads(nThreads);
  for (unsigned iThread = 0; iThread < nThreads; ++iThread) {
    fock_threads[iThread] = std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        fock_threads[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Determine prescreening threshold.
  double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                          this->_lrscf[J]->getSysSettings().basis.integralThreshold);

  if (prescreeningThreshold == 0) {
    double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
    double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
    prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
  }

  double maxDens = 0.0;
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].array().abs().maxCoeff());
    }
  }
  Eigen::MatrixXd maxDensBlocks = this->getShellWiseMaxDens(J, *densityMatrices);
  auto maxDensPtr = maxDensBlocks.data();
  // System I.
  auto basisControllerI = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_I = basisControllerI->getNBasisFunctions();
  // System J.
  auto basisControllerJ = this->_lrscf[J]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_J = basisControllerJ->getNBasisFunctions();
  unsigned nbs_J = basisControllerJ->getReducedNBasisFunctions();
  // aux basis
  auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
  if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
  }
  else if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACCD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
  }
  // Initialize aux basis controller
  auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
  auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);

  // Bases for looper.
  std::shared_ptr<BasisController> superBas = nullptr;
  std::shared_ptr<BasisController> superAuxBas = nullptr;

  if (I == J) {
    superBas = basisControllerI;
    superAuxBas = auxBasisControllerI;
  }
  else {
    BasisFunctionMapper basisMapper(basisControllerI);
    BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
    superBas = basisMapper.getCombinedBasis(basisControllerJ);
    superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
  }
  auto riints = RI_J_IntegralControllerFactory::getInstance().produce(superBas, superAuxBas);
  unsigned nx = superAuxBas->getNBasisFunctions();
  unsigned nxs = superAuxBas->getReducedNBasisFunctions();

  // First contraction.
  std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    sumMat[iSet].resize(this->_nGuess);
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      sumMat[iSet][iGuess] = Eigen::MatrixXd::Zero(nx, nThreads);
    }
  }

  // Distribute function for first contraction.
  auto distribute1 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_J + j;
    unsigned long ji = j * nb_J + i;
    double coul = (i == j ? 0.5 : 1.0) * integral;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto D = (*densityMatrices)[iSet][iGuess].data();
        sumMat[iSet][iGuess](P, threadId) += (D[ij] + D[ji]) * coul;
      }
    }
  };

  // Prescreening function for first contraction.
  auto prescreen1 = [&](unsigned i, unsigned j, unsigned int, double schwarz) {
    unsigned long ij = i * nbs_J + j;
    return (maxDensPtr[ij] * schwarz < prescreeningThreshold);
  };

  // Perform first contraction.
  if (I == J) {
    riints->loopOver3CInts(distribute1, prescreen1);
  }
  else {
    TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisControllerJ, superAuxBas, prescreeningThreshold);
    looper1.loopNoDerivative(distribute1, prescreen1);
  }

  // Solve linear systems.
  std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    coefficients[iSet].resize(this->_nGuess);
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      coefficients[iSet][iGuess] = riints->getLLTMetric().solve(sumMat[iSet][iGuess].rowwise().sum()).eval();
    }
  }

  // Get max values for prescreening.
  Eigen::VectorXd coeffMax = Eigen::VectorXd::Zero(nxs);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto coeffs = superAuxBas->shellWiseAbsMax(coefficients[iSet][iGuess]);
      for (unsigned iShell = 0; iShell < nxs; ++iShell) {
        coeffMax(iShell) = std::max(coeffs(iShell), coeffMax(iShell));
      }
    }
  }

  // Second contraction.
  // Distribute function for second contraction.
  auto distribute2 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_I + j;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = fock_threads[threadId][iSet][iGuess].data();

        // Only perform half the needed contractions and symmetrize afterwards.
        F[ij] += integral * coefficients[iSet][iGuess](P);
      }
    }
  };

  // Prescreening function for second contraction.
  auto prescreen2 = [&](unsigned, unsigned, unsigned K, double schwarz) {
    return (coeffMax(K) * schwarz < prescreeningThreshold);
  };

  // Perform second contraction.
  if (I == J) {
    riints->loopOver3CInts(distribute2, prescreen2);
  }
  else {
    TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisControllerI, superAuxBas, prescreeningThreshold);
    looper2.loopNoDerivative(distribute2, prescreen2);
  }

  // Sum over threads.
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      Eigen::Ref<Eigen::MatrixXd> Fc = fock_threads[0][iSet][iGuess];
      for (unsigned i = 1; i < nThreads; i++) {
        Fc += fock_threads[i][iSet][iGuess];
      }
      Eigen::MatrixXd Fc_sym = Fc + Fc.transpose();
      // The diagonal was fine, so it needs to be halved.
      Fc_sym.diagonal() *= 0.5;
      (*fock)[iSet][iGuess] += Fc_sym;
    }
  }

  Timings::timeTaken("LRSCF -   Sigmavector:   RI-J");
  return fock;
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>
RICoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector:   RI-J");

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Thread safety.
  unsigned nThreads = omp_get_max_threads();
  std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> fock_threads(nThreads);
  for (unsigned iThread = 0; iThread < nThreads; ++iThread) {
    fock_threads[iThread] = std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        fock_threads[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Determine prescreening threshold.
  double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                          this->_lrscf[J]->getSysSettings().basis.integralThreshold);

  if (prescreeningThreshold == 0) {
    double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
    double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
    prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
  }

  double maxDens = 0.0;
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].total().array().abs().maxCoeff());
    }
  }

  Eigen::MatrixXd maxDensBlocks = this->getShellWiseMaxDens(J, *densityMatrices).total();
  auto maxDensPtr = maxDensBlocks.data();
  // System I.
  auto basisControllerI = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_I = basisControllerI->getNBasisFunctions();
  // System J.
  auto basisControllerJ = this->_lrscf[J]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_J = basisControllerJ->getNBasisFunctions();
  unsigned nbs_J = basisControllerJ->getReducedNBasisFunctions();
  // aux basis
  auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
  if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
  }
  else if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACCD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
  }
  // Initialize aux basis controller
  auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
  auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);

  // Bases for looper.
  std::shared_ptr<BasisController> superBas = nullptr;
  std::shared_ptr<BasisController> superAuxBas = nullptr;

  if (I == J) {
    superBas = basisControllerI;
    superAuxBas = auxBasisControllerI;
  }
  else {
    BasisFunctionMapper basisMapper(basisControllerI);
    BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
    superBas = basisMapper.getCombinedBasis(basisControllerJ);
    superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
  }
  auto riints = RI_J_IntegralControllerFactory::getInstance().produce(superBas, superAuxBas);
  unsigned nx = superAuxBas->getNBasisFunctions();
  unsigned nxs = superAuxBas->getReducedNBasisFunctions();

  // First contraction.
  std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    sumMat[iSet].resize(this->_nGuess);
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      sumMat[iSet][iGuess] = Eigen::MatrixXd::Zero(nx, nThreads);
    }
  }

  // Distribute function for first contraction.
  auto distribute1 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_J + j;
    unsigned long ji = j * nb_J + i;
    double coul = (i == j ? 0.5 : 1.0) * integral;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto Da = (*densityMatrices)[iSet][iGuess].alpha.data();
        auto Db = (*densityMatrices)[iSet][iGuess].beta.data();
        sumMat[iSet][iGuess](P, threadId) += (Da[ij] + Da[ji] + Db[ij] + Db[ji]) * coul;
      }
    }
  };

  // Prescreening function for first contraction.
  auto prescreen1 = [&](unsigned i, unsigned j, unsigned int, double schwarz) {
    unsigned long ij = i * nbs_J + j;
    return (maxDensPtr[ij] * schwarz < prescreeningThreshold);
  };

  // Perform first contraction.
  if (I == J) {
    riints->loopOver3CInts(distribute1, prescreen1);
  }
  else {
    TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisControllerJ, superAuxBas, prescreeningThreshold);
    looper1.loopNoDerivative(distribute1, prescreen1);
  }

  // Solve linear systems.
  std::vector<std::vector<Eigen::VectorXd>> coefficients(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    coefficients[iSet].resize(this->_nGuess);
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      coefficients[iSet][iGuess] = riints->getLLTMetric().solve(sumMat[iSet][iGuess].rowwise().sum()).eval();
    }
  }

  // Get max values for prescreening.
  Eigen::VectorXd coeffMax = Eigen::VectorXd::Zero(nxs);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto coeffs = superAuxBas->shellWiseAbsMax(coefficients[iSet][iGuess]);
      for (unsigned iShell = 0; iShell < nxs; ++iShell) {
        coeffMax(iShell) = std::max(coeffs(iShell), coeffMax(iShell));
      }
    }
  }

  // Second contraction.
  // Distribute function for second contraction.
  auto distribute2 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_I + j;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = fock_threads[threadId][iSet][iGuess].data();

        // Only perform half the needed contractions and symmetrize afterwards.
        F[ij] += integral * coefficients[iSet][iGuess](P);
      }
    }
  };

  // Prescreening function for second contraction.
  auto prescreen2 = [&](unsigned, unsigned, unsigned K, double schwarz) {
    return (coeffMax(K) * schwarz < prescreeningThreshold);
  };

  // Perform second contraction.
  if (I == J) {
    riints->loopOver3CInts(distribute2, prescreen2);
  }
  else {
    TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisControllerI, superAuxBas, prescreeningThreshold);
    looper2.loopNoDerivative(distribute2, prescreen2);
  }

  // Sum over threads.
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      Eigen::Ref<Eigen::MatrixXd> Fc = fock_threads[0][iSet][iGuess];
      for (unsigned i = 1; i < nThreads; i++) {
        Fc += fock_threads[i][iSet][iGuess];
      }
      Eigen::MatrixXd Fc_sym = Fc + Fc.transpose();
      // The diagonal was fine, so it needs to be halved.
      Fc_sym.diagonal() *= 0.5;
      (*fock)[iSet][iGuess].alpha += Fc_sym;
      (*fock)[iSet][iGuess].beta += Fc_sym;
    }
  }

  Timings::timeTaken("LRSCF -   Sigmavector:   RI-J");
  return fock;
}

template class RICoulombSigmavector<Options::SCF_MODES::RESTRICTED>;
template class RICoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
