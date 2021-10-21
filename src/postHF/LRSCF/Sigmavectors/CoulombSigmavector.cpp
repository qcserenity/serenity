/**
 * @file CoulombSigmavector.cpp
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
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"

/* Include Serenity Internal Headers */
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
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

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
CoulombSigmavector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
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

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

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
  unsigned ns = maxDensBlocks.rows();

  Timings::takeTime("LRSCF -   Sigmavector:      J");
  if (I == J) {
    // Distribute function II.
    auto distributeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
      unsigned long ij = i * nb_I + j;
      unsigned long ji = j * nb_I + i;
      unsigned long kl = k * nb_I + l;
      unsigned long lk = l * nb_I + k;

      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto D = (*densityMatrices)[iSet][iGuess].data();
          auto F = fock_threads[threadId][iSet][iGuess].data();
          double coul1 = (D[kl] + D[lk]) * integral;
          double coul2 = (D[ij] + D[ji]) * integral;
          F[ij] += coul1;
          F[ji] += coul1;
          F[kl] += coul2;
          F[lk] += coul2;
        }
      }
    };

    // Prescreening function II.
    auto prescreeningFuncII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
      if (maxDens * schwarz < prescreeningThreshold) {
        return true;
      }
      double t1 = maxDensPtr[i * ns + j];
      double t2 = maxDensPtr[j * ns + i];
      double t3 = maxDensPtr[k * ns + l];
      double t4 = maxDensPtr[l * ns + k];
      double maxDBlock = std::max({t1, t2, t3, t4});
      return (maxDBlock * schwarz < prescreeningThreshold);
    };
    TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(), prescreeningThreshold);
    looper.loopNoDerivative(distributeII, prescreeningFuncII, maxDens, nullptr, true);
  }
  else if (I != J) {
    // Distribute function IJ.
    auto distributeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
      unsigned long ij = i * nb_I + j;
      unsigned long ji = j * nb_I + i;
      unsigned long kl = k * nb_J + l;
      unsigned long lk = l * nb_J + k;

      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;

      double coul = perm * integral;
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto D = (*densityMatrices)[iSet][iGuess].data();
          auto F = fock_threads[threadId][iSet][iGuess].data();

          double coul1 = (D[kl] + D[lk]) * coul;
          F[ij] += coul1;
          F[ji] += coul1;
        }
      }
    };

    // Prescreening function IJ.
    auto prescreeningFuncIJ = [&](unsigned, unsigned, unsigned k, unsigned l, double schwarz) {
      if (maxDens * schwarz < prescreeningThreshold) {
        return true;
      }
      double t1 = maxDensPtr[k * ns + l];
      double t2 = maxDensPtr[l * ns + k];
      double maxDBlock = std::max({t1, t2});
      return (maxDBlock * schwarz < prescreeningThreshold);
    };
    CoulombInteractionIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                       this->_lrscf[J]->getBasisController(), prescreeningThreshold);
    looper.loopNoDerivative(distributeIJ, prescreeningFuncIJ);
  }
  Timings::timeTaken("LRSCF -   Sigmavector:      J");

  // Sum over threads.
  for (unsigned iThread = 0; iThread < nThreads; ++iThread) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        (*fock)[iSet][iGuess] += fock_threads[iThread][iSet][iGuess];
      }
    }
  }

  return fock;
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>
CoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
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

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

  double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                          this->_lrscf[J]->getSysSettings().basis.integralThreshold);

  if (prescreeningThreshold == 0) {
    double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
    double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
    prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
  }

  double maxDens = 0.0;
  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].total().array().abs().maxCoeff());
    }
  }

  Eigen::MatrixXd maxDensBlocks = this->getShellWiseMaxDens(J, *densityMatrices).total();
  auto maxDensPtr = maxDensBlocks.data();
  unsigned ns = maxDensBlocks.rows();

  Timings::takeTime("LRSCF -   Sigmavector:      J");
  if (I == J) {
    // Distribute function II.
    auto distributeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
      unsigned long ij = i * nb_I + j;
      unsigned long ji = j * nb_I + i;
      unsigned long kl = k * nb_I + l;
      unsigned long lk = l * nb_I + k;

      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto Da = (*densityMatrices)[iSet][iGuess].alpha.data();
          auto Db = (*densityMatrices)[iSet][iGuess].beta.data();
          auto F = fock_threads[threadId][iSet][iGuess].data();

          double coul1 = (Da[kl] + Da[lk] + Db[kl] + Db[lk]) * integral;
          double coul2 = (Da[ij] + Da[ji] + Db[ij] + Db[ji]) * integral;

          F[ij] += coul1;
          F[ji] += coul1;
          F[kl] += coul2;
          F[lk] += coul2;
        }
      }
    };

    // Prescreening function II.
    auto prescreeningFuncII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
      if (maxDens * schwarz < prescreeningThreshold) {
        return true;
      }
      double t1 = maxDensPtr[i * ns + j];
      double t2 = maxDensPtr[j * ns + i];
      double t3 = maxDensPtr[k * ns + l];
      double t4 = maxDensPtr[l * ns + k];
      double maxDBlock = std::max({t1, t2, t3, t4});
      return (maxDBlock * schwarz < prescreeningThreshold);
    };
    TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(), prescreeningThreshold);
    looper.loopNoDerivative(distributeII, prescreeningFuncII, maxDens, nullptr, true);
  }
  else if (I != J) {
    // Distribute function IJ.
    auto distributeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
      unsigned long ij = i * nb_I + j;
      unsigned long ji = j * nb_I + i;
      unsigned long kl = k * nb_J + l;
      unsigned long lk = l * nb_J + k;

      double perm = 1.0;
      perm *= (i == j) ? 0.5 : 1.0;
      perm *= (k == l) ? 0.5 : 1.0;

      double coul = perm * integral;
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto Da = (*densityMatrices)[iSet][iGuess].alpha.data();
          auto Db = (*densityMatrices)[iSet][iGuess].beta.data();
          auto F = fock_threads[threadId][iSet][iGuess].data();

          double coul1 = (Da[kl] + Da[lk] + Db[kl] + Db[lk]) * coul;
          F[ij] += coul1;
          F[ji] += coul1;
        }
      }
    };

    // Prescreening function IJ.
    auto prescreeningFuncIJ = [&](unsigned, unsigned, unsigned k, unsigned l, double schwarz) {
      if (maxDens * schwarz < prescreeningThreshold) {
        return true;
      }
      double t1 = maxDensPtr[k * ns + l];
      double t2 = maxDensPtr[l * ns + k];
      double maxDBlock = std::max({t1, t2});
      return (maxDBlock * schwarz < prescreeningThreshold);
    };
    CoulombInteractionIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                       this->_lrscf[J]->getBasisController(), prescreeningThreshold);
    looper.loopNoDerivative(distributeIJ, prescreeningFuncIJ);
  }
  Timings::timeTaken("LRSCF -   Sigmavector:      J");

  // Sum over threads.
  for (unsigned iThread = 0; iThread < nThreads; ++iThread) {
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        (*fock)[iSet][iGuess] += fock_threads[iThread][iSet][iGuess];
      }
    }
  }

  return fock;
}

template class CoulombSigmavector<Options::SCF_MODES::RESTRICTED>;
template class CoulombSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
