/**
 * @file DipoleIntegrals.cpp
 * @author Niklas Niemeyer
 *
 * @date Dec. 17, 2018
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
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "grid/GridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DipoleIntegrals<SCFMode>::DipoleIntegrals(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, Point gaugeOrigin)
  : _lrscf(lrscf), _gaugeOrigin(gaugeOrigin) {
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getLengths() {
  if (!_lengths) {
    this->computeIntegrals();
  }
  return _lengths;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getVelocities() {
  if (!_velocities) {
    this->computeIntegrals();
  }
  return _velocities;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<const Eigen::MatrixXd> DipoleIntegrals<SCFMode>::getMagnetics() {
  if (!_magnetics) {
    this->computeIntegrals();
  }
  return _magnetics;
}

template<Options::SCF_MODES SCFMode>
void DipoleIntegrals<SCFMode>::computeIntegrals() {
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> len(_lrscf.size());
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> mag(_lrscf.size());
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> vel(_lrscf.size());

  Timings::takeTime("LRSCF -      Dipole Integrals");
  for (unsigned iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto ints = _lrscf[iSys]->getSys()->getOneElectronIntegralController();
    len[iSys] = ints->getDipoleLengths(_gaugeOrigin);
    vel[iSys] = ints->getDipoleVelocities(_gaugeOrigin);
    mag[iSys] = ints->getDipoleMagnetics(_gaugeOrigin);
  }
  Timings::timeTaken("LRSCF -      Dipole Integrals");

  // Perform AO -> MO transformation and assign shared pointer
  _lengths = std::make_shared<const Eigen::MatrixXd>(ao2mo(len));
  _velocities = std::make_shared<const Eigen::MatrixXd>(ao2mo(vel));
  _magnetics = std::make_shared<const Eigen::MatrixXd>(ao2mo(mag));
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd DipoleIntegrals<SCFMode>::ao2mo(std::vector<std::vector<MatrixInBasis<RESTRICTED>>>& ao_xyz) {
  // Determine length of dipole matrices
  unsigned int nDimension = 0;
  for (unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto nOccupied = _lrscf[iSys]->getNOccupied();
    auto nVirtual = _lrscf[iSys]->getNVirtual();
    for_spin(nOccupied, nVirtual) {
      nDimension += nOccupied_spin * nVirtual_spin;
    };
  }

  Eigen::MatrixXd dipoles = Eigen::MatrixXd::Zero(nDimension, 3);

  unsigned int iStartSys = 0;
  for (unsigned int iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto coeff = _lrscf[iSys]->getCoefficients();
    auto nOccupied = _lrscf[iSys]->getNOccupied();
    auto nVirtual = _lrscf[iSys]->getNVirtual();
    for (unsigned int iMu = 0; iMu < 3; ++iMu) {
      unsigned int iStartSpin = 0;
      for_spin(coeff, nOccupied, nVirtual) {
        Eigen::MatrixXd tmp = coeff_spin.block(0, 0, coeff_spin.rows(), nOccupied_spin).transpose() *
                              ao_xyz[iSys][iMu] * coeff_spin.block(0, nOccupied_spin, coeff_spin.rows(), nVirtual_spin);
        tmp.transposeInPlace();
        dipoles.block(iStartSys + iStartSpin, iMu, nOccupied_spin * nVirtual_spin, 1) =
            Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.cols() * tmp.rows());
        iStartSpin += nOccupied_spin * nVirtual_spin;
      };
    }
    for_spin(nOccupied, nVirtual) {
      iStartSys += nOccupied_spin * nVirtual_spin;
    };
  }

  return dipoles;
}

template class DipoleIntegrals<Options::SCF_MODES::RESTRICTED>;
template class DipoleIntegrals<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
