/**
 * @file DeltaESigmaVector.cpp
 *
 * @date Dec 12, 2018
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
#include "postHF/LRSCF/SigmaVectors/DeltaESigmaVector.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DeltaESigmaVector<SCFMode>::DeltaESigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                              std::vector<Eigen::MatrixXd> b, const double densityScreeningThreshold)
  : SigmaVector<SCFMode>(lrscf, b, densityScreeningThreshold) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
DeltaESigmaVector<SCFMode>::calcF(unsigned int I, unsigned int J,
                                  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  (void)densityMatrices;

  // Orbital-energy differences only contributing if I == J for canonical orbitals
  // LMO-TDDFT in the narrow sense only implemented for supermolecular TDDFT
  if (I == J) {
    // No actual pseudo Fock matrix needed since no AO2MO transformation necessary
    return unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(new std::vector<std::vector<MatrixInBasis<SCFMode>>>());
  }
  else {
    // No contribution here
    return nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
void DeltaESigmaVector<SCFMode>::addToSigma(unsigned int I, unsigned int iStartI,
                                            std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_IJ) {
  Timings::takeTime("LRSCF -    DeltaE Sigmavector");
  // Default pointer to pseudo Fock matrix if contributing
  if (F_IJ) {
    auto orbitalEnergies = this->_lrscf[I]->getEigenvalues();
    auto nOccupied = this->_lrscf[I]->getNOccupied();
    auto nVirtual = this->_lrscf[I]->getNVirtual();

    // Build sigma vector for non-canonical orbitals (LMO-TDDFT)
    if (this->_lrscf[I]->getLRSCFSettings().localMO) {
      // gets the non canonical fock matrix
      auto fockNonCanon = this->_lrscf[I]->getFockNonCanon();

      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
          unsigned int iStartSpin = 0;
          auto& focknon = *fockNonCanon;
          SpinPolarizedData<SCFMode, Eigen::MatrixXd> temp;
          for_spin(nOccupied, nVirtual, focknon, temp) {
            // Sort guess vector in a matrix (occ x virt) such that
            // simple matrix multiplications can be used
            Eigen::MatrixXd guess(nOccupied_spin, nVirtual_spin);
            for (unsigned int i = 0, ia = iStartSpin; i < nOccupied_spin; ++i) {
              for (unsigned int a = 0; a < nVirtual_spin; ++a, ++ia) {
                guess(i, a) = this->_b[iSet][I](ia, iGuess);
              }
            }
            // Perform multiplication of guess vector with off-diagonal Lagrange multiplier
            temp_spin = guess * focknon_spin.bottomRightCorner(nVirtual_spin, nVirtual_spin).transpose();
            temp_spin -= focknon_spin.topLeftCorner(nOccupied_spin, nOccupied_spin) * guess;
            iStartSpin += nOccupied_spin * nVirtual_spin;
          };
          iStartSpin = 0;
          // Go back to vector representation
          for_spin(nOccupied, nVirtual, temp) {
            Eigen::VectorXd guess(nOccupied_spin * nVirtual_spin);
            for (unsigned int i = 0, ia = iStartSpin; i < nOccupied_spin; ++i) {
              for (unsigned int a = 0; a < nVirtual_spin; ++a, ++ia) {
                guess(ia - iStartSpin) = temp_spin(i, a);
              }
            }
            // Set up sigma vector
            this->_sigma[iSet].block(iStartI + iStartSpin, iGuess, nOccupied_spin * nVirtual_spin, 1) += guess;
            iStartSpin += nOccupied_spin * nVirtual_spin;
          };
        }
      }
      // Build sigma vector for canonical orbitals
    }
    else {
      // Prepare vector of orbital-energy differences for each spin
      SpinPolarizedData<SCFMode, Eigen::VectorXd> orbitalEnergyDifferences;
      for_spin(nOccupied, nVirtual, orbitalEnergies, orbitalEnergyDifferences) {
        orbitalEnergyDifferences_spin.resize(nOccupied_spin * nVirtual_spin);
        orbitalEnergyDifferences_spin.setZero();
        for (unsigned int ia = 0; ia < nOccupied_spin * nVirtual_spin; ++ia) {
          unsigned int i = floor(ia / nVirtual_spin);
          unsigned int a = nOccupied_spin + ia - i * nVirtual_spin;
          orbitalEnergyDifferences_spin(ia) = orbitalEnergies_spin(a) - orbitalEnergies_spin(i);
        }
      };

      for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
        unsigned int iStartSpin = 0;
        for_spin(nOccupied, nVirtual, orbitalEnergyDifferences) {
          this->_sigma[iSet].middleRows(iStartI + iStartSpin, nOccupied_spin * nVirtual_spin) +=
              orbitalEnergyDifferences_spin.asDiagonal() *
              this->_b[iSet][I].middleRows(iStartSpin, nOccupied_spin * nVirtual_spin);
          iStartSpin += nOccupied_spin * nVirtual_spin;
        };
      }
      F_IJ.reset();
    }
  }
  else {
    return;
  }
  Timings::timeTaken("LRSCF -    DeltaE Sigmavector");
}

template class DeltaESigmaVector<Options::SCF_MODES::RESTRICTED>;
template class DeltaESigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
