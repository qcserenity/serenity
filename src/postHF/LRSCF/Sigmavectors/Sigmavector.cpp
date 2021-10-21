/**
 * @file Sigmavector.cpp
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers
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
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "basis/Basis.h"
#include "data/OrbitalPair.h"
#include "data/PAOController.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"
/* Include Serenity Internal Headers */
#include "misc/Timing.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Sigmavector<SCFMode>::Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b)
  : _lrscf(lrscf), _nSub(lrscf.size()), _nGuess(b[0].cols()), _nSet(b.size()), _hasBeenCalculated(false) {
  // Construct guess vectors for each subsystem for easier access (b holds guess vectors of total system,
  // i.e. the full subsystem response problem).
  _b.resize(_nSet);
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    _b[iSet].resize(_nSub);
    unsigned int iCount = 0;
    for (unsigned int I = 0; I < _nSub; ++I) {
      // Get position of guess vector for subsystem I
      auto nOccI = _lrscf[I]->getNOccupied();
      auto nVirtI = _lrscf[I]->getNVirtual();
      unsigned int nDimI = 0;
      for_spin(nOccI, nVirtI) {
        nDimI += nOccI_spin * nVirtI_spin;
        iCount += nOccI_spin * nVirtI_spin;
      };
      unsigned int iStartI = iCount - nDimI;
      _b[iSet][I] = b[iSet].block(iStartI, 0, nDimI, _nGuess);
    }
  }
  // Set dimensions for sigma vector and null
  _sigma.resize(b.size());
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    _sigma[iSet].resize(b[iSet].rows(), _nGuess);
    _sigma[iSet].setZero();
  }
}

template<Options::SCF_MODES SCFMode>
Sigmavector<SCFMode>::Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf)
  : _lrscf(lrscf), _nSub(lrscf.size()), _nGuess(1), _nSet(1), _hasBeenCalculated(false) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> Sigmavector<SCFMode>::calcP(unsigned int J) {
  Timings::takeTime("LRSCF -  AO-MO transformation");
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J(
      new std::vector<std::vector<MatrixInBasis<SCFMode>>>(_nSet));
  auto C = _lrscf[J]->getCoefficients();
  auto no = _lrscf[J]->getNOccupied();
  auto nv = _lrscf[J]->getNVirtual();
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      (*P_J)[iSet].emplace_back(_lrscf[J]->getBasisController());
      auto& P = (*P_J)[iSet][iGuess];
      unsigned iStartSpin = 0;
      for_spin(P, C, no, nv) {
        Eigen::Map<Eigen::MatrixXd> b(_b[iSet][J].col(iGuess).data() + iStartSpin, nv_spin, no_spin);
        P_spin = C_spin.middleCols(no_spin, nv_spin) * b * C_spin.middleCols(0, no_spin).transpose();
        P_spin.transposeInPlace();
        iStartSpin += nv_spin * no_spin;
      };
    }
  }
  Timings::timeTaken("LRSCF -  AO-MO transformation");
  return P_J;
}

template<Options::SCF_MODES SCFMode>
void Sigmavector<SCFMode>::addToSigma(unsigned int I, unsigned int iStartI,
                                      std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_I) {
  Timings::takeTime("LRSCF -  AO-MO transformation");
  auto C = _lrscf[I]->getCoefficients();
  auto no = _lrscf[I]->getNOccupied();
  auto nv = _lrscf[I]->getNVirtual();
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      auto& F = (*F_I)[iSet][iGuess];
      unsigned iStartSpin = 0;
      for_spin(F, C, no, nv) {
        F_spin = C_spin.middleCols(no_spin, nv_spin).transpose() * F_spin.transpose() * C_spin.middleCols(0, no_spin);
        _sigma[iSet].block(iStartI + iStartSpin, iGuess, nv_spin * no_spin, 1) +=
            Eigen::Map<Eigen::VectorXd>(F_spin.data(), F_spin.size());
        F_spin.resize(0, 0);
        iStartSpin += nv_spin * no_spin;
      };
    }
  }
  F_I.reset();
  Timings::timeTaken("LRSCF -  AO-MO transformation");
}

template<Options::SCF_MODES SCFMode>
void Sigmavector<SCFMode>::calcSigma() {
  // Set sigma vectors for each set to zero.
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    _sigma[iSet].setZero();
  }
  // Calculate matrix vector product M_IJ*b_J for each set.
  unsigned iCount = 0;
  for (unsigned I = 0; I < _nSub; ++I) {
    auto nOccI = _lrscf[I]->getNOccupied();
    auto nVirtI = _lrscf[I]->getNVirtual();
    unsigned nDimI = 0;
    for_spin(nOccI, nVirtI) {
      nDimI += nOccI_spin * nVirtI_spin;
      iCount += nOccI_spin * nVirtI_spin;
    };
    unsigned iStartI = iCount - nDimI;

    bool exploitSymmetry = (_lrscf.size() > 1 && !_lrscf[I]->getLRSCFSettings().fullFDEc);
    exploitSymmetry = exploitSymmetry && _lrscf[I]->getLRSCFSettings().partialResponseConstruction;
    for (unsigned J = 0; J < _nSub; ++J) {
      if (exploitSymmetry) {
        if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDA) {
          if (I == J || I > J) {
            continue;
          }
        }
        else if (_lrscf[I]->getLRSCFSettings().method == Options::LR_METHOD::TDDFT) {
          if (I > J) {
            continue;
          }
        }
      }

      auto P_J = calcP(J);

      // Simple significance screening.
      double maxDens = 0.0;
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          maxDens = std::max(maxDens, (*P_J)[iSet][iGuess].total().array().abs().maxCoeff());
        }
      }

      // A simple check to see if this block is even populated with data.
      if (maxDens < 1e-8) {
        continue;
      }

      // Detailed significance screening.
      // ToDo: It happens that a guess vector in the set of guess vectors is significant while others
      //      are not. Insignificant guess vectors should be excluded from the set of guess vectors
      //      and their sigma vector should be set to zero.
      auto F_IJ = this->calcF(I, J, std::move(P_J));
      this->addToSigma(I, iStartI, std::move(F_IJ));
    }
  }
  _hasBeenCalculated = true;
}

template<Options::SCF_MODES SCFMode>
std::vector<MatrixInBasis<SCFMode>> Sigmavector<SCFMode>::getPerturbedFockMatrix() {
  // Calculate matrix vector product M_IJ*b_J for each set.
  unsigned int iCount = 0;
  if (_nSub > 1)
    throw SerenityError("More than one subsystem not supported!");

  auto nOccI = _lrscf[0]->getNOccupied();
  auto nVirtI = _lrscf[0]->getNVirtual();
  unsigned int nDimI = 0;
  for_spin(nOccI, nVirtI) {
    nDimI += nOccI_spin * nVirtI_spin;
    iCount += nOccI_spin * nVirtI_spin;
  };
  auto P_J = calcP(0);

  // Simple significance screening.
  double maxDens = 0.0;
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < _nGuess; ++iGuess) {
      auto& P = (*P_J)[iSet][iGuess];
      for_spin(P) {
        maxDens = std::max(maxDens, P_spin.array().abs().maxCoeff());
      };
    }
  }

  // Detailed significance screening.
  // ToDo: It happens, that a guess vector in the set of guess vectors is significant while others
  //      are not. Insignificant guess vectors should be excluded from the set of guess vectors
  //      and their sigma vector should be set to zero.

  auto F_IJ = this->calcF(0, 0, std::move(P_J));

  auto perturbedFockMatrix = (*F_IJ)[0];
  if ((*F_IJ).size() > 1) {
    for (unsigned int iState = 0; iState < perturbedFockMatrix.size(); iState++) {
      perturbedFockMatrix[iState] += (*F_IJ)[1][iState];
    }
  }

  return perturbedFockMatrix;
}

template<>
SPMatrix<Options::SCF_MODES::RESTRICTED> Sigmavector<Options::SCF_MODES::RESTRICTED>::getShellWiseMaxDens(
    unsigned J, std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>& densityMatrices) {
  SPMatrix<Options::SCF_MODES::RESTRICTED> shellWiseMaxDens;
  const auto& basis = this->_lrscf[J]->getBasisController()->getBasis();
  unsigned nShells = basis.size();
  shellWiseMaxDens = SPMatrix<Options::SCF_MODES::RESTRICTED>(Eigen::MatrixXd::Zero(nShells, nShells));
  for (unsigned iShell = 0; iShell < nShells; ++iShell) {
    unsigned nI = basis[iShell]->size();
    unsigned iStart = this->_lrscf[J]->getBasisController()->extendedIndex(iShell);
    for (unsigned jShell = 0; jShell < nShells; ++jShell) {
      unsigned nJ = basis[jShell]->size();
      unsigned jStart = this->_lrscf[J]->getBasisController()->extendedIndex(jShell);
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& densityMatrix = densityMatrices[iSet][iGuess];
          double thisMax = densityMatrix.block(iStart, jStart, nI, nJ).lpNorm<Eigen::Infinity>();
          shellWiseMaxDens(iShell, jShell) = std::max(shellWiseMaxDens(iShell, jShell), thisMax);
        }
      }
    }
  }

  return shellWiseMaxDens;
}

template<>
SPMatrix<Options::SCF_MODES::UNRESTRICTED> Sigmavector<Options::SCF_MODES::UNRESTRICTED>::getShellWiseMaxDens(
    unsigned J, std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>& densityMatrices) {
  SPMatrix<Options::SCF_MODES::UNRESTRICTED> shellWiseMaxDens;
  const auto& basis = this->_lrscf[J]->getBasisController()->getBasis();
  unsigned nShells = basis.size();
  shellWiseMaxDens = SPMatrix<Options::SCF_MODES::UNRESTRICTED>(Eigen::MatrixXd::Zero(nShells, nShells));
  for (unsigned iShell = 0; iShell < nShells; ++iShell) {
    unsigned nI = basis[iShell]->size();
    unsigned iStart = this->_lrscf[J]->getBasisController()->extendedIndex(iShell);
    for (unsigned jShell = 0; jShell < nShells; ++jShell) {
      unsigned nJ = basis[jShell]->size();
      unsigned jStart = this->_lrscf[J]->getBasisController()->extendedIndex(jShell);
      for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          auto& densityMatrix = densityMatrices[iSet][iGuess];
          double alphaMax = densityMatrix.alpha.block(iStart, jStart, nI, nJ).lpNorm<Eigen::Infinity>();
          double betaMax = densityMatrix.beta.block(iStart, jStart, nI, nJ).lpNorm<Eigen::Infinity>();
          shellWiseMaxDens.alpha(iShell, jShell) = std::max(shellWiseMaxDens.alpha(iShell, jShell), alphaMax);
          shellWiseMaxDens.beta(iShell, jShell) = std::max(shellWiseMaxDens.beta(iShell, jShell), betaMax);
        }
      }
    }
  }

  return shellWiseMaxDens;
}

template class Sigmavector<Options::SCF_MODES::RESTRICTED>;
template class Sigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
