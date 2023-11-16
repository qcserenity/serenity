/**
 * @file   CorrespondingOrbitals.cpp
 *
 * @date   Dec 14, 2017
 * @author Anja Massolle
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
#include "analysis/brokenSymmetry/CorrespondingOrbitals.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "integrals/wrappers/Libint.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CorrespondingOrbitals<SCFMode>::CorrespondingOrbitals(std::shared_ptr<SystemController> sys1,
                                                      std::shared_ptr<SystemController> sys2)
  : CorrespondingOrbitals(sys1->getActiveOrbitalController<SCFMode>()->getCoefficients(),
                          sys2->getActiveOrbitalController<SCFMode>()->getCoefficients(),
                          sys1->getNOccupiedOrbitals<SCFMode>(), sys2->getNOccupiedOrbitals<SCFMode>()){};

template<Options::SCF_MODES SCFMode>
CorrespondingOrbitals<SCFMode>::CorrespondingOrbitals(CoefficientMatrix<SCFMode> coeff1, CoefficientMatrix<SCFMode> coeff2,
                                                      SpinPolarizedData<SCFMode, unsigned int> nOcc1,
                                                      SpinPolarizedData<SCFMode, unsigned int> nOcc2)
  : _nOcc1(nOcc1), _nOcc2(nOcc2) {
  printSectionTitle("Calculate Corresponding Orbitals");
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd _SAO =
      libint.compute1eInts(LIBINT_OPERATOR::overlap, coeff1.getBasisController(), coeff2.getBasisController());

  for_spin(coeff1, coeff2, _coeff1, _coeff2, _nOcc1, _nOcc2, _SMO) {
    _coeff1_spin = (coeff1_spin).leftCols(_nOcc1_spin);
    _coeff2_spin = (coeff2_spin).leftCols(_nOcc2_spin);
    _SMO_spin = (_coeff1_spin).transpose() * _SAO * _coeff2_spin;
  };

  for_spin(_U, _V, _SMO) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_SMO_spin, Eigen::ComputeFullU | Eigen::ComputeFullV);
    _U_spin = svd.matrixU();
    _V_spin = svd.matrixV();
  };
};

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, Eigen::MatrixXd>, SpinPolarizedData<SCFMode, Eigen::MatrixXd>>
CorrespondingOrbitals<SCFMode>::getCorrespondingOrbitals() {
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> correspondingOrb1;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> correspondingOrb2;
  for_spin(_coeff1, _coeff2, correspondingOrb1, correspondingOrb2, _U, _V) {
    correspondingOrb1_spin = _coeff1_spin * _U_spin;
    correspondingOrb2_spin = _coeff2_spin * _V_spin;
  };
  return std::make_pair(correspondingOrb1, correspondingOrb2);
};

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd> CorrespondingOrbitals<SCFMode>::getSPrimePrime() {
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> sPrimePrime;
  for_spin(sPrimePrime, _SMO, _U, _V) {
    sPrimePrime_spin = (_U_spin).transpose() * _SMO_spin * _V_spin;
  };
  return sPrimePrime;
};

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, double> CorrespondingOrbitals<SCFMode>::getOverlap(SpinPolarizedData<SCFMode, unsigned int> n) {
  SpinPolarizedData<SCFMode, double> overlap;
  auto sMO = this->getSPrimePrime();
  for_spin(sMO, n, overlap) {
    overlap_spin = (sMO_spin)(n_spin - 1, n_spin - 1);
  };
  return overlap;
};

template class CorrespondingOrbitals<Options::SCF_MODES::RESTRICTED>;
template class CorrespondingOrbitals<Options::SCF_MODES::UNRESTRICTED>;

}; // namespace Serenity
