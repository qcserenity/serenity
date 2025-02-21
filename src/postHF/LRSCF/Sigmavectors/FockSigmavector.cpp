/**
 * @file FockSigmavector.cpp
 *
 * @date Dec 12, 2018
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
#include "postHF/LRSCF/Sigmavectors/FockSigmavector.h"
/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FockSigmavector<SCFMode>::FockSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                          std::vector<Eigen::MatrixXd> b)
  : Sigmavector<SCFMode>(lrscf, b) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
FockSigmavector<SCFMode>::calcF(unsigned I, unsigned J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>) {
  return (I == J) ? std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>() : nullptr;
}

template<Options::SCF_MODES SCFMode>
void FockSigmavector<SCFMode>::addToSigma(unsigned I, unsigned iStartI,
                                          std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_IJ) {
  Timings::takeTime("LRSCF -   Sigmavector:   Fock");

  if (F_IJ) {
    auto no = this->_lrscf[I]->getNOccupied();
    auto nv = this->_lrscf[I]->getNVirtual();
    auto fockPtr = this->_lrscf[I]->getMOFockMatrix();
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        unsigned iStart = 0;
        auto& fock = *fockPtr;
        for_spin(nv, no, fock) {
          Eigen::Map<Eigen::MatrixXd> guess(this->_b[iSet][I].col(iGuess).data() + iStart, nv_spin, no_spin);
          Eigen::Map<Eigen::MatrixXd> sigma(this->_sigma[iSet].col(iGuess).data() + iStartI + iStart, nv_spin, no_spin);
          sigma += fock_spin.bottomRightCorner(nv_spin, nv_spin) * guess - guess * fock_spin.topLeftCorner(no_spin, no_spin);
          iStart += nv_spin * no_spin;
        };
      }
    }
  }

  Timings::timeTaken("LRSCF -   Sigmavector:   Fock");
}

template class FockSigmavector<Options::SCF_MODES::RESTRICTED>;
template class FockSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
