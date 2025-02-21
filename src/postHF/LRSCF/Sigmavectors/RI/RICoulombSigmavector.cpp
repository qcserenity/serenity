/**
 * @file RICoulombSigmavector.cpp
 *
 * @date Dec 07, 2018
 * @author Niklas Niemeyer, Johannes Toelle
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
#include "geometry/Geometry.h"
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

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
RICoulombSigmavector<SCFMode>::calcF(unsigned int I, unsigned int J,
                                     std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector:   RI-J");

  double distance = this->_lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*this->_lrscf[J]->getSys()->getGeometry());
  if (distance >= this->_lrscf[I]->getLRSCFSettings().approxCoulomb[0] && I != J) {
    Timings::timeTaken("LRSCF -   Sigmavector:   RI-J");
    return nullptr;
  }

  // Calculate shell-wise max coefficients in set of perturbed density matrices.
  this->setShellWiseMaxDens(J, (*densityMatrices));
  auto maxDensPtr = this->_maxDensMat.data();

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Auxiliary basis.
  auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
  if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
  }
  else if (this->_lrscf[I]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACCD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
  }

  // System I.
  auto basisControllerI = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
  unsigned nb_I = basisControllerI->getNBasisFunctions();

  // System J.
  auto basisControllerJ = this->_lrscf[J]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);
  unsigned nb_J = basisControllerJ->getNBasisFunctions();
  unsigned nbs_J = basisControllerJ->getReducedNBasisFunctions();

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

  // Thread safety.
  std::vector<std::vector<std::vector<MatrixInBasis<RESTRICTED>>>> Fc(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    Fc[iThread] = std::vector<std::vector<MatrixInBasis<RESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Fc[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Use symmetry.
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> D_sym(
      this->_nSet, std::vector<MatrixInBasis<RESTRICTED>>(this->_nGuess, this->_lrscf[J]->getBasisController()));
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto D = (*densityMatrices)[iSet][iGuess].total();
      D_sym[iSet][iGuess] += D;
      D_sym[iSet][iGuess] += D.transpose();
    }
  }

  // First contraction.
  std::vector<std::vector<Eigen::MatrixXd>> sumMat(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    sumMat[iSet].resize(this->_nGuess);
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      sumMat[iSet][iGuess] = Eigen::MatrixXd::Zero(nx, this->_nThreads);
    }
  }

  auto distribute1 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_J + j;
    double coul = (i == j ? 0.5 : 1.0) * integral;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto D = D_sym[iSet][iGuess].data();
        sumMat[iSet][iGuess](P, threadId) += D[ij] * coul;
      }
    }
  };

  auto prescreen1 = [&](unsigned i, unsigned j, unsigned int, double schwarz) {
    unsigned long ij = i * nbs_J + j;
    return (maxDensPtr[ij] * schwarz < this->_prescreeningThreshold);
  };

  // Perform first contraction.
  TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisControllerJ, superAuxBas, this->_prescreeningThreshold);
  looper1.loopNoDerivative(distribute1, prescreen1);

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
  auto distribute2 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
    unsigned long ij = i * nb_I + j;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = Fc[threadId][iSet][iGuess].data();
        F[ij] += integral * coefficients[iSet][iGuess](P);
      }
    }
  };

  auto prescreen2 = [&](unsigned, unsigned, unsigned P, double schwarz) {
    return (coeffMax(P) * schwarz < this->_prescreeningThreshold);
  };

  // Perform second contraction.
  TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisControllerI, superAuxBas, this->_prescreeningThreshold);
  looper2.loopNoDerivative(distribute2, prescreen2);
  Timings::timeTaken("LRSCF -   Sigmavector:   RI-J");

  // Sum over threads.
  Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      auto& F = (*fock)[iSet][iGuess];
      for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
        F += Fc[iThread][iSet][iGuess];
      }
      for_spin(F) {
        F_spin += F_spin.transpose().eval();
        F_spin.diagonal() *= 0.5;
      };
    }
  }
  Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");

  return fock;
}

template class RICoulombSigmavector<Options::SCF_MODES::RESTRICTED>;
template class RICoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
