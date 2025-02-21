/**
 * @file CoulombSigmavector.cpp
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
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CoulombSigmavector<SCFMode>::CoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                std::vector<Eigen::MatrixXd> b)
  : Sigmavector<SCFMode>(lrscf, b) {
}

template<Options::SCF_MODES SCFMode>
CoulombSigmavector<SCFMode>::CoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : Sigmavector<SCFMode>(lrscf) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
CoulombSigmavector<SCFMode>::calcF(unsigned int I, unsigned int J,
                                   std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector:      J");

  double distance = this->_lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*this->_lrscf[J]->getSys()->getGeometry());
  if (distance >= this->_lrscf[I]->getLRSCFSettings().approxCoulomb[0] && I != J) {
    Timings::timeTaken("LRSCF -   Sigmavector:      J");
    return nullptr;
  }

  // Calculate shell-wise max coefficients in set of perturbed density matrices.
  this->setShellWiseMaxDens(J, (*densityMatrices));

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

  // Prescreening.
  auto maxDensPtr = this->_maxDensMat.data();
  unsigned ns = this->_maxDensMat.rows();

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

  // Distribute function II.
  auto distributeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ij = i * nb_I + j;
    unsigned kl = k * nb_I + l;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = Fc[iThread][iSet][iGuess].data();
        auto D = D_sym[iSet][iGuess].data();
        F[ij] += D[kl] * integral;
        F[kl] += D[ij] * integral;
      }
    }
  };

  // Distribute function IJ.
  auto distributeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ij = i * nb_I + j;
    unsigned kl = k * nb_J + l;

    double perm = 1.0;
    perm *= (i == j) ? 0.5 : 1.0;
    perm *= (k == l) ? 0.5 : 1.0;

    double coul = perm * integral;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = Fc[iThread][iSet][iGuess].data();
        auto D = D_sym[iSet][iGuess].data();
        F[ij] += D[kl] * coul;
      }
    }
  };

  // Prescreening function II.
  auto prescreeningFuncII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    if (this->_maxDens * schwarz < this->_prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + j];
    double t2 = maxDensPtr[j * ns + i];
    double t3 = maxDensPtr[k * ns + l];
    double t4 = maxDensPtr[l * ns + k];
    double maxDBlock = std::max({t1, t2, t3, t4});
    return (maxDBlock * schwarz < this->_prescreeningThreshold);
  };

  // Prescreening function IJ.
  auto prescreeningFuncIJ = [&](unsigned, unsigned, unsigned k, unsigned l, double schwarz) {
    if (this->_maxDens * schwarz < this->_prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[k * ns + l];
    double t2 = maxDensPtr[l * ns + k];
    double maxDBlock = std::max({t1, t2});
    return (maxDBlock * schwarz < this->_prescreeningThreshold);
  };

  if (I == J) {
    TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                      this->_prescreeningThreshold);
    looper.loopNoDerivative(distributeII, prescreeningFuncII, this->_maxDens, nullptr, true);
  }
  else if (I != J) {
    CoulombInteractionIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                       this->_lrscf[J]->getBasisController(), this->_prescreeningThreshold);
    looper.loopNoDerivative(distributeIJ, prescreeningFuncIJ);
  }
  Timings::timeTaken("LRSCF -   Sigmavector:      J");

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
      };
    }
  }
  Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");

  return fock;
}

template class CoulombSigmavector<Options::SCF_MODES::RESTRICTED>;
template class CoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
