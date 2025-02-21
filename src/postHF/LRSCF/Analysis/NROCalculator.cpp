/**
 * @file NROCalculator.cpp
 *
 * @date Apr 23, 2021
 * @author Anton Rikus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "postHF/LRSCF/Analysis/NROCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
NROCalculator<SCFMode>::NROCalculator(std::vector<Eigen::MatrixXd> solutionvectors,
                                      std::shared_ptr<LRSCFController<SCFMode>> lrscfController)
  : _lrscfcontroller(lrscfController), _XY(solutionvectors){};

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> NROCalculator<SCFMode>::getNROs(unsigned int iFreq) {
  auto nocc = _lrscfcontroller->getNOccupied();
  auto nvirt = _lrscfcontroller->getNVirtual();
  auto& coeffmat = _lrscfcontroller->getCoefficients();
  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> NROs;
  Eigen::MatrixXd xpy = _XY[0].middleCols(iFreq * 3, 3);
  for (unsigned j = 0; j < 3; j++) {
    unsigned ia_start = 0;
    for_spin(coeffmat, nocc, nvirt, NROs) {
      Eigen::MatrixXd xpymatrix =
          (Eigen::Map<Eigen::MatrixXd>(xpy.col(j).data() + ia_start, nvirt_spin, nocc_spin)).transpose();
      Eigen::JacobiSVD<Eigen::MatrixXd> svdp(xpymatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::MatrixXd U = (svdp.matrixU());
      Eigen::MatrixXd V = (svdp.matrixV());
      for_spin(_singularValues, nocc) {
        _singularValues_spin.push_back(Eigen::MatrixXd::Zero(nocc_spin, 3));
        _singularValues_spin[iFreq].col(j) = svdp.singularValues() / svdp.singularValues().sum();
      };
      NROs_spin.push_back(coeffmat_spin.leftCols(nocc_spin) * U);
      NROs_spin.push_back(coeffmat_spin.middleCols(nocc_spin, nvirt_spin) * V.leftCols(nocc_spin));
      ia_start += nocc_spin * nvirt_spin;
    };
  }
  return NROs;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> NROCalculator<SCFMode>::getSingularValues() {
  return _singularValues;
}

template class NROCalculator<Options::SCF_MODES::RESTRICTED>;
template class NROCalculator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */