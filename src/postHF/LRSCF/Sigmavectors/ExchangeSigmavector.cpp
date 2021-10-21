/**
 * @file ExchangeSigmavector.cpp
 *
 * @date Dec 11, 2018
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
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"

/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "integrals/looper/ExchangeInteractionIntLooper.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "io/FormattedOutputStream.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/EmbeddingOptions.h"
#include "settings/Settings.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ExchangeSigmavector<SCFMode>::ExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                  std::vector<Eigen::MatrixXd> b, const std::vector<int> pm,
                                                  bool densFitK, bool densFitLRK)
  : Sigmavector<SCFMode>(lrscf, b), _pm(pm), _densFitK(densFitK), _densFitLRK(densFitLRK) {
}

template<Options::SCF_MODES SCFMode>
ExchangeSigmavector<SCFMode>::ExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : Sigmavector<SCFMode>(lrscf), _pm({0}) {
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
ExchangeSigmavector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
  this->setParameters(I, J, false);

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
  unsigned ns = maxDensBlocks.rows();

  // Exchange ratio for lambda functions.
  double exchangeRatio = 0.0;

  // Distribute function II.
  auto distributeExchangeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nb_I + k;
    unsigned long il = i * nb_I + l;
    unsigned long jl = j * nb_I + l;
    unsigned long jk = j * nb_I + k;
    unsigned long ki = k * nb_I + i;
    unsigned long kj = k * nb_I + j;
    unsigned long li = l * nb_I + i;
    unsigned long lj = l * nb_I + j;
    double exc = integral * exchangeRatio;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      double exc_B = _pm[iSet] * exc;
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = fock_threads[threadId][iSet][iGuess].data();
        auto D = (*densityMatrices)[iSet][iGuess].data();

        // Contractions from A exchange term.
        F[ik] += D[jl] * exc;
        F[il] += D[jk] * exc;
        F[jk] += D[il] * exc;
        F[jl] += D[ik] * exc;
        F[ki] += D[lj] * exc;
        F[kj] += D[li] * exc;
        F[li] += D[kj] * exc;
        F[lj] += D[ki] * exc;

        // Contractions from B exchange term.
        if (_pm[iSet]) {
          F[ik] += D[lj] * exc_B;
          F[il] += D[kj] * exc_B;
          F[jk] += D[li] * exc_B;
          F[jl] += D[ki] * exc_B;
          F[ki] += D[jl] * exc_B;
          F[kj] += D[il] * exc_B;
          F[li] += D[jk] * exc_B;
          F[lj] += D[ik] * exc_B;
        }
      }
    }
  };

  // Prescreening function II.
  auto prescreeningFuncII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    if (maxDens * schwarz < prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double t5 = maxDensPtr[k * ns + i];
    double t6 = maxDensPtr[k * ns + j];
    double t7 = maxDensPtr[l * ns + i];
    double t8 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (maxDBlock * schwarz < prescreeningThreshold);
  };

  // Distribute function IJ.
  auto distributeExchangeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nb_I + k;
    unsigned long jl = j * nb_J + l;
    unsigned long lj = l * nb_J + j;

    double exc = integral * exchangeRatio;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      double exc_B = _pm[iSet] * exc;
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto D = (*densityMatrices)[iSet][iGuess].data();
        auto F = fock_threads[threadId][iSet][iGuess].data();

        // Contractions from A exchange term.
        F[ik] += D[jl] * exc;

        // Contractions from B exchange term.
        if (_pm[iSet]) {
          F[ik] += D[lj] * exc_B;
        }
      }
    }
  };

  // Prescreening function IJ.
  auto prescreeningFuncIJ = [&](unsigned, unsigned j, unsigned, unsigned l, double schwarz) {
    if (maxDens * schwarz < prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[j * ns + l];
    double t2 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2});
    return (maxDBlock * schwarz < prescreeningThreshold);
  };

  // Loop and contract integrals (conventional exchange)!
  if (_hfExchangeRatio > 0.0 && !_densFitK) {
    OutputControl::dOut << " **** Performing No RI-K ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:      K");
    exchangeRatio = _hfExchangeRatio;
    if (I == J) {
      TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold);
      looper.loopNoDerivative(distributeExchangeII, prescreeningFuncII, maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), prescreeningThreshold);
      looper.loopNoDerivative(distributeExchangeIJ, prescreeningFuncIJ);
    }
    Timings::timeTaken("LRSCF -   Sigmavector:      K");
  }

  // Loop and contract integrals (long-range exchange)!
  if (_lrExchangeRatio > 0.0 && !_densFitLRK) {
    OutputControl::dOut << " **** Performing No RI-ErfK ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:    LRK");
    exchangeRatio = _lrExchangeRatio;
    if (I == J) {
      TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::erf_coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold, _mu);
      looper.loopNoDerivative(distributeExchangeII, prescreeningFuncII, maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(LIBINT_OPERATOR::erf_coulomb, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), prescreeningThreshold, _mu);
      looper.loopNoDerivative(distributeExchangeIJ, prescreeningFuncIJ);
    }
    Timings::timeTaken("LRSCF -   Sigmavector:    LRK");
  }

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
ExchangeSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
  this->setParameters(I, J, false);

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Thread safety.
  unsigned nThreads = omp_get_max_threads();
  std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> fock_threads(nThreads);
  for (unsigned iThread = 0; iThread < nThreads; ++iThread) {
    fock_threads[iThread] = std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        fock_threads[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

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
  unsigned ns = maxDensBlocks.rows();

  // Exchange ratio for lambda functions.
  double exchangeRatio = 0.0;

  // Distribute function II.
  auto distributeExchangeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nb_I + k;
    unsigned long il = i * nb_I + l;
    unsigned long jl = j * nb_I + l;
    unsigned long jk = j * nb_I + k;
    unsigned long ki = k * nb_I + i;
    unsigned long kj = k * nb_I + j;
    unsigned long li = l * nb_I + i;
    unsigned long lj = l * nb_I + j;
    double exc = integral * exchangeRatio;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      double exc_B = _pm[iSet] * exc;
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto Fa = fock_threads[threadId][iSet][iGuess].alpha.data();
        auto Fb = fock_threads[threadId][iSet][iGuess].beta.data();
        auto Da = (*densityMatrices)[iSet][iGuess].alpha.data();
        auto Db = (*densityMatrices)[iSet][iGuess].beta.data();

        // Contractions from A exchange term.
        Fa[ik] += Da[jl] * exc;
        Fa[il] += Da[jk] * exc;
        Fa[jk] += Da[il] * exc;
        Fa[jl] += Da[ik] * exc;
        Fa[ki] += Da[lj] * exc;
        Fa[kj] += Da[li] * exc;
        Fa[li] += Da[kj] * exc;
        Fa[lj] += Da[ki] * exc;

        Fb[ik] += Db[jl] * exc;
        Fb[il] += Db[jk] * exc;
        Fb[jk] += Db[il] * exc;
        Fb[jl] += Db[ik] * exc;
        Fb[ki] += Db[lj] * exc;
        Fb[kj] += Db[li] * exc;
        Fb[li] += Db[kj] * exc;
        Fb[lj] += Db[ki] * exc;

        // Contractions from B exchange term.
        if (_pm[iSet]) {
          Fa[ik] += Da[lj] * exc_B;
          Fa[il] += Da[kj] * exc_B;
          Fa[jk] += Da[li] * exc_B;
          Fa[jl] += Da[ki] * exc_B;
          Fa[ki] += Da[jl] * exc_B;
          Fa[kj] += Da[il] * exc_B;
          Fa[li] += Da[jk] * exc_B;
          Fa[lj] += Da[ik] * exc_B;

          Fb[ik] += Db[lj] * exc_B;
          Fb[il] += Db[kj] * exc_B;
          Fb[jk] += Db[li] * exc_B;
          Fb[jl] += Db[ki] * exc_B;
          Fb[ki] += Db[jl] * exc_B;
          Fb[kj] += Db[il] * exc_B;
          Fb[li] += Db[jk] * exc_B;
          Fb[lj] += Db[ik] * exc_B;
        }
      }
    }
  };

  // Prescreening function II.
  auto prescreeningFuncII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    if (maxDens * schwarz < prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[i * ns + k];
    double t2 = maxDensPtr[i * ns + l];
    double t3 = maxDensPtr[j * ns + k];
    double t4 = maxDensPtr[j * ns + l];
    double t5 = maxDensPtr[k * ns + i];
    double t6 = maxDensPtr[k * ns + j];
    double t7 = maxDensPtr[l * ns + i];
    double t8 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (maxDBlock * schwarz < prescreeningThreshold);
  };

  // Distribute function IJ.
  auto distributeExchangeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned threadId) {
    unsigned long ik = i * nb_I + k;
    unsigned long jl = j * nb_J + l;
    unsigned long lj = l * nb_J + j;

    double exc = integral * exchangeRatio;
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      double exc_B = exc * _pm[iSet];
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto Da = (*densityMatrices)[iSet][iGuess].alpha.data();
        auto Db = (*densityMatrices)[iSet][iGuess].beta.data();
        auto Fa = fock_threads[threadId][iSet][iGuess].alpha.data();
        auto Fb = fock_threads[threadId][iSet][iGuess].beta.data();

        // Contractions from A exchange term.
        Fa[ik] += Da[jl] * exc;
        Fb[ik] += Db[jl] * exc;

        // Contractions from B exchange term.
        if (_pm[iSet]) {
          Fa[ik] += Da[lj] * exc_B;
          Fb[ik] += Db[lj] * exc_B;
        }
      }
    }
  };

  // Prescreening function IJ.
  auto prescreeningFuncIJ = [&](unsigned, unsigned j, unsigned, unsigned l, double schwarz) {
    if (maxDens * schwarz < prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[j * ns + l];
    double t2 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2});
    return (maxDBlock * schwarz < prescreeningThreshold);
  };

  // Loop and contract integrals (conventional exchange)!
  if (_hfExchangeRatio > 0.0 && !_densFitK) {
    OutputControl::dOut << " **** Performing No RI-K ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:      K");
    exchangeRatio = _hfExchangeRatio;
    if (I == J) {
      TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold);
      looper.loopNoDerivative(distributeExchangeII, prescreeningFuncII, maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(LIBINT_OPERATOR::coulomb, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), prescreeningThreshold);
      looper.loopNoDerivative(distributeExchangeIJ, prescreeningFuncIJ);
    }
    Timings::timeTaken("LRSCF -   Sigmavector:      K");
  }

  // Loop and contract integrals (long-range exchange)!
  if (_lrExchangeRatio > 0.0 && !_densFitLRK) {
    OutputControl::dOut << " **** Performing No RI-ErfK ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:    LRK");
    exchangeRatio = _lrExchangeRatio;
    if (I == J) {
      TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::erf_coulomb, 0, this->_lrscf[I]->getBasisController(),
                                        prescreeningThreshold, _mu);
      looper.loopNoDerivative(distributeExchangeII, prescreeningFuncII, maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(LIBINT_OPERATOR::erf_coulomb, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), prescreeningThreshold, _mu);
      looper.loopNoDerivative(distributeExchangeIJ, prescreeningFuncIJ);
    }
    Timings::timeTaken("LRSCF -   Sigmavector:    LRK");
  }

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

template<Options::SCF_MODES SCFMode>
void ExchangeSigmavector<SCFMode>::setParameters(unsigned I, unsigned J, bool rpaScreen) {
  if (I == J && this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    Functional funcI = resolveFunctional(this->_lrscf[I]->getLRSCFSettings().func);
    if (this->_lrscf[I]->getLRSCFSettings().func == CompositeFunctionals::XCFUNCTIONALS::NONE) {
      funcI = resolveFunctional(this->_lrscf[I]->getSysSettings().dft.functional);
    }
    if (!rpaScreen) {
      _hfExchangeRatio = funcI.getHfExchangeRatio();
      _lrExchangeRatio = funcI.getLRExchangeRatio();
      _mu = funcI.getRangeSeparationParameter();
    }
  }
  else if (I != J) {
    auto funcNadd = resolveFunctional(this->_lrscf[I]->getLRSCFSettings().embedding.naddXCFunc);
    _hfExchangeRatio = funcNadd.getHfExchangeRatio();
    _lrExchangeRatio = funcNadd.getLRExchangeRatio();
    _mu = funcNadd.getRangeSeparationParameter();
    auto embList = this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList;
    if (embList.size() > 0) {
      if (embList.size() <= I || embList.size() <= J) {
        throw SerenityError("embeddingModeList.size() smaller than the number of subsystems");
      }
      // Look if exact or approximate embedding.
      if ((embList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT && embList[J] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
          (embList[I] == Options::KIN_EMBEDDING_MODES::HOFFMANN && embList[J] == Options::KIN_EMBEDDING_MODES::HOFFMANN) ||
          (embList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA && embList[J] == Options::KIN_EMBEDDING_MODES::HUZINAGA) ||
          (embList[I] == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA &&
           embList[J] == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA)) {
        auto funcNaddExact = resolveFunctional(this->_lrscf[I]->getLRSCFSettings().embedding.naddXCFuncList[0]);
        _hfExchangeRatio = funcNaddExact.getHfExchangeRatio();
        _lrExchangeRatio = funcNaddExact.getLRExchangeRatio();
        _mu = funcNaddExact.getRangeSeparationParameter();
      }
      else {
        auto funcNaddApprox = resolveFunctional(this->_lrscf[I]->getLRSCFSettings().embedding.naddXCFuncList[1]);
        _hfExchangeRatio = funcNaddApprox.getHfExchangeRatio();
        _lrExchangeRatio = funcNaddApprox.getLRExchangeRatio();
        _mu = funcNaddApprox.getRangeSeparationParameter();
      }
    }
  }
}

template class ExchangeSigmavector<Options::SCF_MODES::RESTRICTED>;
template class ExchangeSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
