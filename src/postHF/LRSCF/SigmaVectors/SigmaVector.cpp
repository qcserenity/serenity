/**
 * @file SigmaVector.cpp
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
#include "postHF/LRSCF/SigmaVectors/SigmaVector.h"
#include "data/OrbitalPair.h"
#include "data/PAOController.h"
#include "data/SingleSubstitution.h"
/* Include Serenity Internal Headers */
#include "misc/Timing.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SigmaVector<SCFMode>::SigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                  std::vector<Eigen::MatrixXd> b, const double densityScreeningThreshold)
  : _lrscf(lrscf),
    _densityScreeningThreshold(densityScreeningThreshold),
    _nSub(lrscf.size()),
    _nGuess(b[0].cols()),
    _nSet(b.size()),
    _hasBeenCalculated(false) {
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
SigmaVector<SCFMode>::SigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                  const double densityScreeningThreshold)
  : _lrscf(lrscf),
    _densityScreeningThreshold(densityScreeningThreshold),
    _nSub(lrscf.size()),
    _nGuess(1),
    _nSet(1),
    _hasBeenCalculated(false) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> SigmaVector<SCFMode>::calcP(unsigned int I) {
  Timings::takeTime("LRSCF -        Pseudo density");
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_I(
      new std::vector<std::vector<MatrixInBasis<SCFMode>>>(_nSet));
  auto coeff = _lrscf[I]->getCoefficients();
  auto nOccupied = _lrscf[I]->getNOccupied();
  auto nVirtual = _lrscf[I]->getNVirtual();
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < _nGuess; ++iGuess) {
      (*P_I)[iSet].emplace_back(_lrscf[I]->getBasisController());
      auto& P = (*P_I)[iSet][iGuess];
      unsigned int iStartSpin = 0;
      for_spin(P, coeff, nOccupied, nVirtual) {
        Eigen::MatrixXd guess(nOccupied_spin, nVirtual_spin);
        for (unsigned int i = 0, ia = iStartSpin; i < nOccupied_spin; ++i) {
          for (unsigned int a = 0; a < nVirtual_spin; ++a, ++ia) {
            guess(i, a) = _b[iSet][I](ia, iGuess);
          }
        }
        P_spin =
            coeff_spin.leftCols(nOccupied_spin) * guess * coeff_spin.middleCols(nOccupied_spin, nVirtual_spin).transpose();
        iStartSpin += nOccupied_spin * nVirtual_spin;
      };
    }
  }
  Timings::timeTaken("LRSCF -        Pseudo density");
  return P_I;
}

template<Options::SCF_MODES SCFMode>
void SigmaVector<SCFMode>::addToSigma(unsigned int I, unsigned int iStartI,
                                      std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_IJ) {
  Timings::takeTime("LRSCF -  AO2MO transformation");
  auto coeff = _lrscf[I]->getCoefficients();
  auto nOccupied = _lrscf[I]->getNOccupied();
  auto nVirtual = _lrscf[I]->getNVirtual();
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < _nGuess; ++iGuess) {
      auto& pF = (*F_IJ)[iSet][iGuess];
      unsigned int iStartSpin = 0;
      for_spin(pF, coeff, nOccupied, nVirtual) {
        pF_spin = coeff_spin.leftCols(nOccupied_spin).transpose() * pF_spin *
                  coeff_spin.middleCols(nOccupied_spin, nVirtual_spin);
        pF_spin.transposeInPlace();
        _sigma[iSet].block(iStartI + iStartSpin, iGuess, nOccupied_spin * nVirtual_spin, 1) +=
            Eigen::Map<Eigen::VectorXd>(pF_spin.data(), pF_spin.cols() * pF_spin.rows());
        pF_spin.resize(0, 0);
        iStartSpin += nOccupied_spin * nVirtual_spin;
      };
    }
  }
  F_IJ.reset();
  Timings::timeTaken("LRSCF -  AO2MO transformation");
}

template<Options::SCF_MODES SCFMode>
void SigmaVector<SCFMode>::calcSigma() {
  // Set sigma vectors for each set to zero
  for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
    _sigma[iSet].setZero();
  }
  // Calculate matrix vector product M_IJ*b_J for each set
  unsigned int iCount = 0;
  for (unsigned int I = 0; I < _nSub; ++I) {
    auto nOccI = _lrscf[I]->getNOccupied();
    auto nVirtI = _lrscf[I]->getNVirtual();
    unsigned int nDimI = 0;
    for_spin(nOccI, nVirtI) {
      nDimI += nOccI_spin * nVirtI_spin;
      iCount += nOccI_spin * nVirtI_spin;
    };
    unsigned int iStartI = iCount - nDimI;
    for (unsigned int J = 0; J < _nSub; ++J) {
      auto P_J = calcP(J);

      // simple significance screening
      double maxDens = 0.0;
      for (unsigned int iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < _nGuess; ++iGuess) {
          auto& P = (*P_J)[iSet][iGuess];
          for_spin(P) {
            maxDens = std::max(maxDens, P_spin.array().abs().maxCoeff());
          };
        }
      }
      // If negligible: Skip the contributions arising from guess vectors belonging to this subsystem
      if (maxDens < _densityScreeningThreshold)
        continue;

      // detailed significance screening
      // ToDo: It happens that a guess vector in the set of guess vectors is significant while others
      //      are not. Insignificant guess vectors should be excluded from the set of guess vectors
      //      and their sigma vector should be set to zero.
      auto F_IJ = this->calcF(I, J, std::move(P_J));
      this->addToSigma(I, iStartI, std::move(F_IJ));
    }
  }
  _hasBeenCalculated = true;
}

template class SigmaVector<Options::SCF_MODES::RESTRICTED>;
template class SigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
