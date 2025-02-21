/**
 * @file MRICoulombSigmavector.cpp
 *
 * @date Nov 01, 2022
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/Sigmavectors/RI/MRICoulombSigmavector.h"
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
MRICoulombSigmavector<SCFMode>::MRICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                      std::vector<Eigen::MatrixXd> b)
  : Sigmavector<SCFMode>(lrscf, b) {
}

template<Options::SCF_MODES SCFMode>
MRICoulombSigmavector<SCFMode>::MRICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : Sigmavector<SCFMode>(lrscf) {
}

template<Options::SCF_MODES SCFMode>
void MRICoulombSigmavector<SCFMode>::calcSigma() {
  Timings::takeTime("LRSCF -   Sigmavector:  MRI-J");

  // Auxiliary basis.
  auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
  if (this->_lrscf[0]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
  }
  else if (this->_lrscf[0]->getLRSCFSettings().densFitJ == Options::DENS_FITS::ACCD) {
    basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
  }

  // Determine number of basis functions.
  long nxsuper = 0;
  Eigen::VectorXi nxstart = Eigen::VectorXi::Zero(this->_nSub);
  for (unsigned I = 0; I < this->_nSub; ++I) {
    unsigned nx = this->_lrscf[I]->getBasisController(basispurpose)->getNBasisFunctions();
    nxsuper += nx;
    if (I < this->_nSub - 1) {
      nxstart[I + 1] = nxstart[I] + nx;
    }
  }

  std::vector<Eigen::MatrixXd> X(this->_nSet, Eigen::MatrixXd::Zero(nxsuper, this->_nGuess));
  for (unsigned I = 0; I < this->_nSub; ++I) {
    auto basisController = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
    auto auxBasisController = this->_lrscf[I]->getBasisController(basispurpose);
    auto riints = RI_J_IntegralControllerFactory::getInstance().produce(basisController, auxBasisController);

    unsigned nb = basisController->getNBasisFunctions();
    unsigned nbs = basisController->getReducedNBasisFunctions();
    unsigned nx = auxBasisController->getNBasisFunctions();
    double prescreeningThreshold = basisController->getPrescreeningThreshold();

    auto P = this->calcP(I);
    std::vector<std::vector<MatrixInBasis<RESTRICTED>>> D_sym(
        this->_nSet, std::vector<MatrixInBasis<RESTRICTED>>(this->_nGuess, basisController));
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto D = (*P)[iSet][iGuess].total();
        D_sym[iSet][iGuess] += D;
        D_sym[iSet][iGuess] += D.transpose();
      }
    }

    std::vector<std::vector<Eigen::MatrixXd>> sumMat(
        this->_nSet, std::vector<Eigen::MatrixXd>(this->_nGuess, Eigen::MatrixXd::Zero(nx, this->_nThreads)));

    auto distribute1 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
      unsigned long ij = i * nb + j;
      double coul = (i == j ? 0.5 : 1.0) * integral;
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto D = D_sym[iSet][iGuess].data();
          sumMat[iSet][iGuess](P, threadId) += D[ij] * coul;
        }
      }
    };

    this->setShellWiseMaxDens(I, (*P));
    auto maxDensPtr = this->_maxDensMat.data();
    auto prescreen1 = [&](unsigned i, unsigned j, unsigned int, double schwarz) {
      unsigned long ij = i * nbs + j;
      return (maxDensPtr[ij] * schwarz < prescreeningThreshold);
    };

    TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisController, auxBasisController, prescreeningThreshold);
    looper1.loopNoDerivative(distribute1, prescreen1);

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Eigen::VectorXd rhs = sumMat[iSet][iGuess].rowwise().sum();
        X[iSet].col(iGuess).segment(nxstart[I], nx) = riints->getLLTMetric().solve(rhs);
      }
    }
  }

  auto& libint = Libint::getInstance();
  std::vector<Eigen::MatrixXd> Z(this->_nSet, Eigen::MatrixXd::Zero(nxsuper, this->_nGuess));
  for (unsigned I = 0; I < this->_nSub; ++I) {
    auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
    unsigned nxI = auxBasisControllerI->getNBasisFunctions();

    for (unsigned J = I + 1; J < this->_nSub; ++J) {
      auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);
      unsigned nxJ = auxBasisControllerJ->getNBasisFunctions();

      if (this->skipThisInteraction(I, J) || I == J) {
        continue;
      }

      double distance =
          this->_lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*this->_lrscf[J]->getSys()->getGeometry());
      if (distance < this->_lrscf[I]->getLRSCFSettings().approxCoulomb[0] ||
          distance >= this->_lrscf[I]->getLRSCFSettings().approxCoulomb[1]) {
        continue;
      }

      Eigen::MatrixXd S = libint.compute1eInts(LIBINT_OPERATOR::coulomb, auxBasisControllerJ, auxBasisControllerI);
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        Z[iSet].middleRows(nxstart[I], nxI).noalias() += S * X[iSet].middleRows(nxstart[J], nxJ);
        Z[iSet].middleRows(nxstart[J], nxJ).noalias() += S.transpose() * X[iSet].middleRows(nxstart[I], nxI);
      }
    }
  }

  long iaStart = 0;
  for (unsigned I = 0; I < this->_nSub; ++I) {
    auto no = this->_lrscf[I]->getNOccupied();
    auto nv = this->_lrscf[I]->getNVirtual();
    auto C = this->_lrscf[I]->getCoefficients();

    auto basisController = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
    auto auxBasisController = this->_lrscf[I]->getBasisController(basispurpose);
    auto riints = RI_J_IntegralControllerFactory::getInstance().produce(basisController, auxBasisController);

    unsigned nb = basisController->getNBasisFunctions();
    unsigned nx = auxBasisController->getNBasisFunctions();
    unsigned nxs = auxBasisController->getReducedNBasisFunctions();
    double prescreeningThreshold = basisController->getPrescreeningThreshold();

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      Eigen::Ref<Eigen::MatrixXd> rhs = Z[iSet].middleRows(nxstart[I], nx);
      X[iSet].middleRows(nxstart[I], nx) = riints->getLLTMetric().solve(rhs);
    }

    std::vector<std::vector<std::vector<Eigen::MatrixXd>>> Fc(this->_nThreads);
    for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
      Fc[iThread] = std::vector<std::vector<Eigen::MatrixXd>>(this->_nSet);
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          Fc[iThread][iSet].push_back(Eigen::MatrixXd::Zero(nb, nb));
        }
      }
    }

    // Second contraction.
    auto distribute2 = [&](unsigned i, unsigned j, unsigned P, double integral, unsigned threadId) {
      unsigned long ij = i * nb + j;
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto F = Fc[threadId][iSet][iGuess].data();
          F[ij] += integral * X[iSet](nxstart[I] + P, iGuess);
        }
      }
    };

    Eigen::VectorXd coeffMax = Eigen::VectorXd::Zero(nxs);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto coeffs = auxBasisController->shellWiseAbsMax(X[iSet].middleRows(nxstart[I], nx).col(iGuess));
        for (unsigned iShell = 0; iShell < nxs; ++iShell) {
          coeffMax(iShell) = std::max(coeffs(iShell), coeffMax(iShell));
        }
      }
    }

    auto prescreen2 = [&](unsigned, unsigned, unsigned P, double schwarz) {
      return (coeffMax(P) * schwarz < prescreeningThreshold);
    };

    TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisController, auxBasisController, prescreeningThreshold);
    looper2.loopNoDerivative(distribute2, prescreen2);

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        // Sum over threads.
        for (unsigned iThread = 1; iThread < this->_nThreads; ++iThread) {
          Fc[0][iSet][iGuess] += Fc[iThread][iSet][iGuess];
        }

        Eigen::MatrixXd F = Fc[0][iSet][iGuess] + Fc[0][iSet][iGuess].transpose();
        F.diagonal() *= 0.5;

        // AO -> MO transformation.
        unsigned iStartSpin = 0;
        for_spin(C, no, nv) {
          F = (C_spin.middleCols(no_spin, nv_spin).transpose() * F.transpose() * C_spin.middleCols(0, no_spin)).eval();
          this->_sigma[iSet].block(iaStart + iStartSpin, iGuess, nv_spin * no_spin, 1) +=
              Eigen::Map<Eigen::VectorXd>(F.data(), F.size());
          iStartSpin += nv_spin * no_spin;
        };
      }
    }

    for_spin(no, nv) {
      iaStart += nv_spin * no_spin;
    };
  }

  this->_hasBeenCalculated = true;
  Timings::timeTaken("LRSCF -   Sigmavector:  MRI-J");
}

template class MRICoulombSigmavector<Options::SCF_MODES::RESTRICTED>;
template class MRICoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
