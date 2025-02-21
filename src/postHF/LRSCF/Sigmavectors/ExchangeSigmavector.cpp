/**
 * @file ExchangeSigmavector.cpp
 *
 * @date Dec 11, 2018
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
ExchangeSigmavector<SCFMode>::ExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                  const std::vector<int> pm, bool densFitK, bool densFitLRK)
  : Sigmavector<SCFMode>(lrscf), _pm(pm), _densFitK(densFitK), _densFitLRK(densFitLRK) {
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
ExchangeSigmavector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
  this->setParameters(I, J, false);

  // Calculate shell-wise max coefficients in set of perturbed density matrices.
  this->setShellWiseMaxDens(J, (*densityMatrices));

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Thread safety.
  std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> Fx(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    Fx[iThread] = std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Fx[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

  // Prescreening.
  auto maxDensPtr = this->_maxDensMat.data();
  unsigned ns = this->_maxDensMat.rows();

  // Exchange ratio for lambda functions.
  double exchangeRatio = 0.0;

  // Use symmetry.
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> D_sym(
      this->_nSet, std::vector<MatrixInBasis<RESTRICTED>>(this->_nGuess, this->_lrscf[J]->getBasisController()));
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      D_sym[iSet][iGuess] = (*densityMatrices)[iSet][iGuess];
      D_sym[iSet][iGuess] += _pm[iSet] * (*densityMatrices)[iSet][iGuess].transpose();
    }
  }

  // Distribute function II.
  auto distributeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ii = i * nb_I;
    unsigned jj = j * nb_I;
    unsigned kk = k * nb_I;
    unsigned ll = l * nb_I;
    unsigned ik = ii + k;
    unsigned il = ii + l;
    unsigned jk = jj + k;
    unsigned jl = jj + l;
    unsigned ki = kk + i;
    unsigned kj = kk + j;
    unsigned li = ll + i;
    unsigned lj = ll + j;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = Fx[iThread][iSet][iGuess].data();
        auto D = D_sym[iSet][iGuess].data();
        F[ik] += D[jl] * integral;
        F[il] += D[jk] * integral;
        F[jk] += D[il] * integral;
        F[jl] += D[ik] * integral;
        if (!_pm[iSet]) {
          F[ki] += D[lj] * integral;
          F[li] += D[kj] * integral;
          F[kj] += D[li] * integral;
          F[lj] += D[ki] * integral;
        }
      }
    }
  };

  // Distribute function IJ.
  auto distributeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ik = i * nb_I + k;
    unsigned jl = j * nb_J + l;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto F = Fx[iThread][iSet][iGuess].data();
        auto D = D_sym[iSet][iGuess].data();
        F[ik] += D[jl] * integral;
      }
    }
  };

  // Prescreening function II.
  auto prescreenII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = exchangeRatio * schwarz;
    if (xschwarz * this->_maxDens < this->_prescreeningThreshold) {
      return true;
    }
    unsigned ii = i * ns;
    unsigned jj = j * ns;
    unsigned kk = k * ns;
    unsigned ll = l * ns;
    double t1 = maxDensPtr[ii + k];
    double t2 = maxDensPtr[ii + l];
    double t3 = maxDensPtr[jj + k];
    double t4 = maxDensPtr[jj + l];
    double t5 = maxDensPtr[kk + i];
    double t6 = maxDensPtr[kk + j];
    double t7 = maxDensPtr[ll + i];
    double t8 = maxDensPtr[ll + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (xschwarz * maxDBlock < this->_prescreeningThreshold);
  };

  // Prescreening function IJ.
  auto prescreenIJ = [&](unsigned, unsigned j, unsigned, unsigned l, double schwarz) {
    double xschwarz = exchangeRatio * schwarz;
    if (xschwarz * this->_maxDens < this->_prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[j * ns + l];
    double t2 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2});
    return (xschwarz * maxDBlock < this->_prescreeningThreshold);
  };

  auto contractExchange = [&](double ax, LIBINT_OPERATOR op, double mu) {
    exchangeRatio = ax;

    // Loop and contract.
    if (I == J) {
      TwoElecFourCenterIntLooper looper(op, 0, this->_lrscf[I]->getBasisController(), this->_prescreeningThreshold, mu);
      looper.loopNoDerivative(distributeII, prescreenII, this->_maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(op, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), this->_prescreeningThreshold, mu);
      looper.loopNoDerivative(distributeIJ, prescreenIJ);
    }

    // Sum up fock matrices.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
          (*fock)[iSet][iGuess] += exchangeRatio * Fx[iThread][iSet][iGuess];
          // Set zero if needed another time.
          if (op == LIBINT_OPERATOR::coulomb && _lrExchangeRatio > 0.0) {
            Fx[iThread][iSet][iGuess].setZero();
          }
        }
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  };

  // HF exchange.
  if (_hfExchangeRatio > 0.0 && !_densFitK) {
    Timings::takeTime("LRSCF -   Sigmavector:      K");
    contractExchange(_hfExchangeRatio, LIBINT_OPERATOR::coulomb, 0.0);
    Timings::timeTaken("LRSCF -   Sigmavector:      K");
  }

  // LR exchange.
  if (_lrExchangeRatio > 0.0 && !_densFitLRK) {
    Timings::takeTime("LRSCF -   Sigmavector:    LRK");
    contractExchange(_lrExchangeRatio, LIBINT_OPERATOR::erf_coulomb, _mu);
    Timings::timeTaken("LRSCF -   Sigmavector:    LRK");
  }

  // Symmetrize.
  if (I == J) {
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        (*fock)[iSet][iGuess] += _pm[iSet] * (*fock)[iSet][iGuess].transpose().eval();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  }

  return fock;
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>
ExchangeSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
  this->setParameters(I, J, false);

  // Calculate shell-wise max coefficients in set of perturbed density matrices.
  this->setShellWiseMaxDens(J, (*densityMatrices));

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Thread safety.
  std::vector<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> Fx(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    Fx[iThread] = std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>(this->_nSet);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Fx[iThread][iSet].emplace_back(this->_lrscf[I]->getBasisController());
      }
    }
  }

  // Number of AOs in subsystem I and J.
  unsigned nb_I = this->_lrscf[I]->getBasisController()->getNBasisFunctions();
  unsigned nb_J = this->_lrscf[J]->getBasisController()->getNBasisFunctions();

  // Prescreening.
  auto maxDensPtr = this->_maxDensMat.data();
  unsigned ns = this->_maxDensMat.rows();

  // Exchange ratio for lambda functions.
  double exchangeRatio = 0.0;

  // Use symmetry.
  std::vector<std::vector<MatrixInBasis<UNRESTRICTED>>> D_sym(
      this->_nSet, std::vector<MatrixInBasis<UNRESTRICTED>>(this->_nGuess, this->_lrscf[J]->getBasisController()));
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      D_sym[iSet][iGuess].alpha = (*densityMatrices)[iSet][iGuess].alpha;
      D_sym[iSet][iGuess].beta = (*densityMatrices)[iSet][iGuess].beta;
      D_sym[iSet][iGuess].alpha += _pm[iSet] * (*densityMatrices)[iSet][iGuess].alpha.transpose();
      D_sym[iSet][iGuess].beta += _pm[iSet] * (*densityMatrices)[iSet][iGuess].beta.transpose();
    }
  }

  // Distribute function II.
  auto distributeII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ii = i * nb_I;
    unsigned jj = j * nb_I;
    unsigned kk = k * nb_I;
    unsigned ll = l * nb_I;
    unsigned ik = ii + k;
    unsigned il = ii + l;
    unsigned jk = jj + k;
    unsigned jl = jj + l;
    unsigned ki = kk + i;
    unsigned kj = kk + j;
    unsigned li = ll + i;
    unsigned lj = ll + j;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto Fa = Fx[iThread][iSet][iGuess].alpha.data();
        auto Da = D_sym[iSet][iGuess].alpha.data();
        auto Fb = Fx[iThread][iSet][iGuess].beta.data();
        auto Db = D_sym[iSet][iGuess].beta.data();
        Fa[ik] += Da[jl] * integral;
        Fa[il] += Da[jk] * integral;
        Fa[jk] += Da[il] * integral;
        Fa[jl] += Da[ik] * integral;
        Fb[ik] += Db[jl] * integral;
        Fb[il] += Db[jk] * integral;
        Fb[jk] += Db[il] * integral;
        Fb[jl] += Db[ik] * integral;
        if (!_pm[iSet]) {
          Fa[ki] += Da[lj] * integral;
          Fa[li] += Da[kj] * integral;
          Fa[kj] += Da[li] * integral;
          Fa[lj] += Da[ki] * integral;
          Fb[ki] += Db[lj] * integral;
          Fb[li] += Db[kj] * integral;
          Fb[kj] += Db[li] * integral;
          Fb[lj] += Db[ki] * integral;
        }
      }
    }
  };

  // Distribute function IJ.
  auto distributeIJ = [&](unsigned i, unsigned j, unsigned k, unsigned l, double integral, unsigned iThread) {
    unsigned ik = i * nb_I + k;
    unsigned jl = j * nb_J + l;

    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto Fa = Fx[iThread][iSet][iGuess].alpha.data();
        auto Da = D_sym[iSet][iGuess].alpha.data();
        auto Fb = Fx[iThread][iSet][iGuess].beta.data();
        auto Db = D_sym[iSet][iGuess].beta.data();
        Fa[ik] += Da[jl] * integral;
        Fb[ik] += Db[jl] * integral;
      }
    }
  };

  // Prescreening function II.
  auto prescreenII = [&](unsigned i, unsigned j, unsigned k, unsigned l, double schwarz) {
    double xschwarz = exchangeRatio * schwarz;
    if (xschwarz * this->_maxDens < this->_prescreeningThreshold) {
      return true;
    }
    unsigned ii = i * ns;
    unsigned jj = j * ns;
    unsigned kk = k * ns;
    unsigned ll = l * ns;
    double t1 = maxDensPtr[ii + k];
    double t2 = maxDensPtr[ii + l];
    double t3 = maxDensPtr[jj + k];
    double t4 = maxDensPtr[jj + l];
    double t5 = maxDensPtr[kk + i];
    double t6 = maxDensPtr[kk + j];
    double t7 = maxDensPtr[ll + i];
    double t8 = maxDensPtr[ll + j];
    double maxDBlock = std::max({t1, t2, t3, t4, t5, t6, t7, t8});
    return (xschwarz * maxDBlock < this->_prescreeningThreshold);
  };

  // Prescreening function IJ.
  auto prescreenIJ = [&](unsigned, unsigned j, unsigned, unsigned l, double schwarz) {
    double xschwarz = exchangeRatio * schwarz;
    if (xschwarz * this->_maxDens < this->_prescreeningThreshold) {
      return true;
    }
    double t1 = maxDensPtr[j * ns + l];
    double t2 = maxDensPtr[l * ns + j];
    double maxDBlock = std::max({t1, t2});
    return (xschwarz * maxDBlock < this->_prescreeningThreshold);
  };

  auto contractExchange = [&](double ax, LIBINT_OPERATOR op, double mu) {
    exchangeRatio = ax;

    // Loop and contract.
    if (I == J) {
      TwoElecFourCenterIntLooper looper(op, 0, this->_lrscf[I]->getBasisController(), this->_prescreeningThreshold, mu);
      looper.loopNoDerivative(distributeII, prescreenII, this->_maxDens, nullptr, true);
    }
    else if (I != J) {
      ExchangeInteractionIntLooper looper(op, 0, this->_lrscf[I]->getBasisController(),
                                          this->_lrscf[J]->getBasisController(), this->_prescreeningThreshold, mu);
      looper.loopNoDerivative(distributeIJ, prescreenIJ);
    }

    // Sum up fock matrices.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
          (*fock)[iSet][iGuess].alpha += exchangeRatio * Fx[iThread][iSet][iGuess].alpha;
          (*fock)[iSet][iGuess].beta += exchangeRatio * Fx[iThread][iSet][iGuess].beta;
          // Set zero if needed another time.
          if (op == LIBINT_OPERATOR::coulomb && _lrExchangeRatio > 0.0) {
            Fx[iThread][iSet][iGuess].alpha.setZero();
            Fx[iThread][iSet][iGuess].beta.setZero();
          }
        }
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  };

  // HF exchange.
  if (_hfExchangeRatio > 0.0 && !_densFitK) {
    Timings::takeTime("LRSCF -   Sigmavector:      K");
    contractExchange(_hfExchangeRatio, LIBINT_OPERATOR::coulomb, 0.0);
    Timings::timeTaken("LRSCF -   Sigmavector:      K");
  }

  // LR exchange.
  if (_lrExchangeRatio > 0.0 && !_densFitLRK) {
    Timings::takeTime("LRSCF -   Sigmavector:    LRK");
    contractExchange(_lrExchangeRatio, LIBINT_OPERATOR::erf_coulomb, _mu);
    Timings::timeTaken("LRSCF -   Sigmavector:    LRK");
  }

  // Symmetrize.
  if (I == J) {
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        (*fock)[iSet][iGuess].alpha += _pm[iSet] * (*fock)[iSet][iGuess].alpha.transpose().eval();
        (*fock)[iSet][iGuess].beta += _pm[iSet] * (*fock)[iSet][iGuess].beta.transpose().eval();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  }

  return fock;
}

template<Options::SCF_MODES SCFMode>
void ExchangeSigmavector<SCFMode>::setParameters(unsigned I, unsigned J, bool rpaScreen) {
  if (I == J) {
    if (this->_lrscf[I]->getSysSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      Functional funcI = this->_lrscf[I]->getLRSCFSettings().customFunc.basicFunctionals.size()
                             ? Functional(this->_lrscf[I]->getLRSCFSettings().customFunc)
                             : resolveFunctional(this->_lrscf[I]->getLRSCFSettings().func);
      if (this->_lrscf[I]->getLRSCFSettings().func == CompositeFunctionals::XCFUNCTIONALS::NONE) {
        funcI = this->_lrscf[I]->getSysSettings().customFunc.basicFunctionals.size()
                    ? Functional(this->_lrscf[I]->getSysSettings().customFunc)
                    : resolveFunctional(this->_lrscf[I]->getSysSettings().dft.functional);
      }
      if (!rpaScreen) {
        _hfExchangeRatio = funcI.getHfExchangeRatio();
        _lrExchangeRatio = funcI.getLRExchangeRatio();
        _mu = funcI.getRangeSeparationParameter();
      }
    }
    else {
      _hfExchangeRatio = 1.0;
      _lrExchangeRatio = 0.0;
      _mu = 0.0;
    }
  }
  else {
    auto funcNadd = this->_lrscf[I]->getLRSCFSettings().embedding.customNaddXCFunc.basicFunctionals.size()
                        ? Functional(this->_lrscf[I]->getLRSCFSettings().embedding.customNaddXCFunc)
                        : resolveFunctional(this->_lrscf[I]->getLRSCFSettings().embedding.naddXCFunc);
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
