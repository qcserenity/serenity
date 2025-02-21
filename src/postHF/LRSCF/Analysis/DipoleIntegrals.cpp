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
#include "integrals/OneElectronIntegralController.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/RICC2/CC2Controller.h"
#include "settings/LRSCFOptions.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

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

  long iStart = 0;
  Timings::takeTime("LRSCF -      Dipole Integrals");
  for (unsigned iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto ints = _lrscf[iSys]->getSys()->getOneElectronIntegralController();
    auto nb = _lrscf[iSys]->getBasisController()->getNBasisFunctions();

    len[iSys] = ints->getDipoleLengths(_gaugeOrigin);
    vel[iSys] = ints->getDipoleVelocities(_gaugeOrigin);
    mag[iSys] = ints->getDipoleMagnetics(_gaugeOrigin);

    iStart += nb * nb;
  }
  Timings::timeTaken("LRSCF -      Dipole Integrals");

  // Perform AO -> MO transformation and assign shared pointer
  _lengths = std::make_shared<const Eigen::MatrixXd>(this->ao2mo(len));
  _velocities = std::make_shared<const Eigen::MatrixXd>(this->ao2mo(vel));
  _magnetics = std::make_shared<const Eigen::MatrixXd>(this->ao2mo(mag));
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd DipoleIntegrals<SCFMode>::ao2mo(std::vector<std::vector<MatrixInBasis<RESTRICTED>>>& ao_xyz) {
  Eigen::MatrixXd dipoles(0, 3);

  long iStartSysSpin = 0;
  for (unsigned iSys = 0; iSys < _lrscf.size(); ++iSys) {
    auto no = _lrscf[iSys]->getNOccupied();
    auto nv = _lrscf[iSys]->getNVirtual();
    auto P = _lrscf[iSys]->getParticleCoefficients();
    auto H = _lrscf[iSys]->getHoleCoefficients();
    auto cc2 = _lrscf[iSys]->getCC2Controller();

    auto settings = _lrscf[iSys]->getLRSCFSettings();
    bool pq_space = ((cc2 && (!settings.frequencies.empty() || (settings.ccexdens || settings.cctrdens)) &&
                      settings.method != Options::LR_METHOD::CISD) ||
                     this->_fullSpace);

    for_spin(nv, no, P, H) {
      unsigned nb = nv_spin + no_spin;
      unsigned np = pq_space ? nb : nv_spin;
      unsigned nh = pq_space ? nb : no_spin;
      dipoles.conservativeResize(dipoles.rows() + np * nh, 3);
      for (unsigned j = 0; j < 3; ++j) {
        Eigen::MatrixXd tmp =
            P_spin.middleCols(pq_space ? 0 : no_spin, np).transpose() * ao_xyz[iSys][j] * H_spin.middleCols(0, nh);
        dipoles.block(iStartSysSpin, j, np * nh, 1) = Eigen::Map<Eigen::VectorXd>(tmp.data(), np * nh);
      }
      iStartSysSpin += np * nh;
    };
  }

  return dipoles;
}

template<Options::SCF_MODES SCFMode>
Point DipoleIntegrals<SCFMode>::getGaugeOrigin() {
  return _gaugeOrigin;
}

template class DipoleIntegrals<Options::SCF_MODES::RESTRICTED>;
template class DipoleIntegrals<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
