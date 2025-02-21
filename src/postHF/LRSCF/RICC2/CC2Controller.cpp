/**
 * @file CC2Controller.cpp
 *
 * @date Apr 20, 2020
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/RICC2/CC2Controller.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "math/LaplaceMinimaxWrapper.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CC2Controller<SCFMode>::CC2Controller(std::shared_ptr<LRSCFController<SCFMode>> lrscf)
  : _lrscf(lrscf),
    _hfEnergy(lrscf->getSys()->template getElectronicStructure<SCFMode>()->getEnergy()),
    _e(lrscf->getEigenvalues()),
    _C(lrscf->getCoefficients()),
    _riTrafo(*lrscf->getRIIntegrals()->getAuxTrafoPtr()),
    _no(lrscf->getNOccupied()),
    _nv(lrscf->getNVirtual()),
    _nb(lrscf->getRIIntegrals()->getNBasisFunctions()),
    _nx(lrscf->getRIIntegrals()->getNTransformedAuxBasisFunctions()),
    _nxb(lrscf->getRIIntegrals()->getNAuxBasisFunctions()),
    _sss(lrscf->getLRSCFSettings().sss),
    _oss(lrscf->getLRSCFSettings().oss),
    _soss(_sss + _oss),
    _Jia(lrscf->getRIIntegrals()->getJiaPtr()),
    _Jij(lrscf->getRIIntegrals()->getJijPtr()),
    _integrals(lrscf->getRIIntegrals()->getIntegralPtr()),
    _xwfModel(lrscf->getLRSCFSettings().method),
    _settings(lrscf->getLRSCFSettings()),
    _prescreeningThreshold(lrscf->getBasisController()->getPrescreeningThreshold()) {
  if (_settings.ltconv != 0.0 || !_settings.frequencies.empty()) {
    // Use Laplace Transformation to avoid N5 steps.
    double LUMO = +std::numeric_limits<double>::infinity();
    double HOMO = -std::numeric_limits<double>::infinity();
    double HUMO = -std::numeric_limits<double>::infinity();
    double LOMO = +std::numeric_limits<double>::infinity();
    for_spin(_e, _no, _nv) {
      LUMO = std::min(_e_spin(_no_spin), LUMO);
      HOMO = std::max(_e_spin(_no_spin - 1), HOMO);
      HUMO = std::max(_e_spin(_no_spin + _nv_spin - 1), HUMO);
      LOMO = std::min(_e_spin(0), LOMO);
    };
    getMinimaxRoots(_roots, _weights, 2 * (LUMO - HOMO), 2 * (HUMO - LOMO),
                    _settings.ltconv == 0 ? _settings.conv : _settings.ltconv);
  }
}

template<>
void CC2Controller<RESTRICTED>::initialize() {
  if (_xwfModel == Options::LR_METHOD::CC2) {
    _Bij = std::make_shared<SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>>(_no * _no, _nx);
    _Bia = std::make_shared<SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>>(_nv * _no, _nx);
  }
  else {
    _Bij = _Jij;
    _Bia = _Jia;
  }

  // set helper ints
  _nv_a = _nv;
  _no_a = _no;

  _alpha = _nv_a * _no_a;
  _nDim = _nv_a * _no_a;

  _Wia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
  _Xia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
  _Yia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
  _Zia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);

  _Fai = Eigen::VectorXd::Zero(_no * _nv);
  _Fia = Eigen::VectorXd::Zero(_no * _nv);
  _singles = Eigen::VectorXd::Zero(_no * _nv);
  _residual = Eigen::VectorXd::Zero(_no * _nv);
  _gsLagrange = Eigen::VectorXd::Zero(_no * _nv);

  _gsDensity = Eigen::MatrixXd::Zero(_no + _nv, _no + _nv);

  _eVirt = Eigen::MatrixXd::Zero(_nv, _nv);
  _eVirt.colwise() += _e.segment(_no, _nv);
  _eVirt.rowwise() += _e.segment(_no, _nv).transpose();

  _eOcc = Eigen::MatrixXd::Zero(_no, _no);
  _eOcc.colwise() += _e.segment(0, _no);
  _eOcc.rowwise() += _e.segment(0, _no).transpose();

  this->calculateGroundstate();
  this->calculateE();

  // In the triplet case, set spin parts of antisymmetrized amplitudes.
  if (_settings.triplet) {
    _soss = _sss - _oss;
  }
} /* this->initialize() restricted */

template<>
void CC2Controller<UNRESTRICTED>::initialize() {
  if (_xwfModel == Options::LR_METHOD::CC2) {
    _Bij = std::make_shared<SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>>();
    _Bia = std::make_shared<SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>>();
    auto& Bij = *_Bij;
    auto& Bia = *_Bia;
    for_spin(Bij, Bia, _nv, _no) {
      Bij_spin = Eigen::MatrixXd::Zero(_no_spin * _no_spin, _nx);
      Bia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    };
  }
  else {
    _Bij = _Jij;
    _Bia = _Jia;
  }

  // set helper ints
  _nv_a = _nv.alpha;
  _nv_b = _nv.beta;
  _no_a = _no.alpha;
  _no_b = _no.beta;

  _alpha = _nv_a * _no_a;
  _beta = _nv_b * _no_b;

  _nDim = _alpha + _beta;

  _Wia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  _Xia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  _Yia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  _Zia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  for_spin(_Wia, _Xia, _Yia, _Zia, _nv, _no) {
    _Wia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    _Xia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    _Yia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    _Zia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
  };

  _Fai = Eigen::VectorXd::Zero(_nDim);
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _singles = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);
  _gsLagrange = Eigen::VectorXd::Zero(_nDim);

  for_spin(_gsDensity, _no, _nv) {
    _gsDensity_spin = Eigen::MatrixXd::Zero(_no_spin + _nv_spin, _no_spin + _nv_spin);
  };

  for_spin(_eVirt, _eOcc, _nv, _no, _e) {
    _eVirt_spin = Eigen::MatrixXd::Zero(_nv_spin, _nv_spin);
    _eVirt_spin.colwise() += _e_spin.segment(_no_spin, _nv_spin);
    _eVirt_spin.rowwise() += _e_spin.segment(_no_spin, _nv_spin).transpose();

    _eOcc_spin = Eigen::MatrixXd::Zero(_no_spin, _no_spin);
    _eOcc_spin.colwise() += _e_spin.segment(0, _no_spin);
    _eOcc_spin.rowwise() += _e_spin.segment(0, _no_spin).transpose();
  };

  _eVirtab = Eigen::MatrixXd::Zero(_nv_a, _nv_b);
  _eVirtab.colwise() += _e.alpha.segment(_no_a, _nv_a);
  _eVirtab.rowwise() += _e.beta.segment(_no_b, _nv_b).transpose();

  _eOccab = Eigen::MatrixXd::Zero(_no_a, _no_b);
  _eOccab.colwise() += _e.alpha.segment(0, _no_a);
  _eOccab.rowwise() += _e.beta.segment(0, _no_b).transpose();

  this->calculateGroundstate();
  this->calculateE();
} /* this->initialize() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity