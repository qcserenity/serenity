/**
 * @file XWFController.cpp
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
#include "postHF/LRSCF/Sigmavectors/RICC2/XWFController.h"

/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "math/LaplaceMinimaxWrapper.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
XWFController<SCFMode>::XWFController(std::shared_ptr<LRSCFController<SCFMode>> lrscf)
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
Eigen::MatrixXd XWFController<RESTRICTED>::getAmplitudes(unsigned i, unsigned j, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getAmplitudes(unsigned i, unsigned j, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Bia->alpha.middleRows(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j)) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Bia->beta.middleRows(i * _nv_b, _nv_b) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j)) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j)) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getAmplitudesA(unsigned i, unsigned j, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getAmplitudesA() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getAmplitudesA(unsigned i, unsigned j, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Bia->alpha.middleRows(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j)) - _eVirt.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Bia->beta.middleRows(i * _nv_b, _nv_b) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j)) - _eVirt.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j)) - _eVirtab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getAmplitudesA() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getAmplitudesV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getAmplitudesV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getAmplitudesV(unsigned a, unsigned b, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(a * _no_a, _no_a) * _Bia->alpha.middleRows(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b)) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Bia->beta.middleRows(a * _no_b, _no_b) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b)) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(a * _no_a, _no_a) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b)) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getAmplitudesAV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return amp.cwiseQuotient(denom);
} /* this->getAmplitudesAV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getAmplitudesAV(unsigned a, unsigned b, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(a * _no_a, _no_a) * _Bia->alpha.middleRows(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b)) + _eOcc.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Bia->beta.middleRows(a * _no_b, _no_b) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b)) + _eOcc.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Bia->alpha.middleRows(a * _no_a, _no_a) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b)) + _eOccab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getAmplitudesAV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getRightAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose() +
                        _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getRightAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getRightAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Bia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.alpha.middleRows(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) + eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(i * _nv_b, _nv_b) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Bia->beta.middleRows(i * _nv_b, _nv_b) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) + eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) + eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getRightAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getRightAmplitudesA(unsigned i, unsigned j, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose() +
                        _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getRightAmplitudesA() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getRightAmplitudesA(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Bia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.alpha.middleRows(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) + eigenvalue) - _eVirt.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(i * _nv_b, _nv_b) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Bia->beta.middleRows(i * _nv_b, _nv_b) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) + eigenvalue) - _eVirt.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Bia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Bia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) + eigenvalue) - _eVirtab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getRightAmplitudesA() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getRightAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose() +
                        _Bia->middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getRightAmplitudesV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getRightAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Bia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Bia->alpha.middleRows(a * _no_a, _no_a) * _Xia.alpha.middleRows(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) + eigenvalue) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(a * _no_b, _no_b) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Bia->beta.middleRows(a * _no_b, _no_b) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) + eigenvalue) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Bia->alpha.middleRows(a * _no_a, _no_a) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) + eigenvalue) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getRightAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getRightAmplitudesAV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose() +
                        _Bia->middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return amp.cwiseQuotient(denom);
} /* this->getRightAmplitudesAV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getRightAmplitudesAV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Bia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Bia->alpha.middleRows(a * _no_a, _no_a) * _Xia.alpha.middleRows(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) + eigenvalue) + _eOcc.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(a * _no_b, _no_b) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Bia->beta.middleRows(a * _no_b, _no_b) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) + eigenvalue) + _eOcc.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Bia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Bia->alpha.middleRows(a * _no_a, _no_a) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) + eigenvalue) + _eOccab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getRightAmplitudesAV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getLeftAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Zia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Zia.middleRows(j * _nv, _nv).transpose() +
                        _trafoVector.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _trafoVector.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getLeftAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getLeftAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _trafoVector.segment(i * _nv_a, _nv_a) * _Fia.segment(j * _nv_a, _nv_a).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _trafoVector.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) + eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _trafoVector.segment(_alpha + i * _nv_b, _nv_b) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fia.segment(_alpha + i * _nv_b, _nv_b) * _trafoVector.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) + eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _trafoVector.segment(i * _nv_a, _nv_a) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _trafoVector.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) + eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getLeftAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getLeftAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Zia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Zia.middleRows(b * _no, _no).transpose() +
                        _trafoVector.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _trafoVector.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getLeftAmplitudesV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getLeftAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _trafoVector.segment(a * _no_a, _no_a) * _Fia.segment(b * _no_a, _no_a).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _trafoVector.segment(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) + eigenvalue) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _trafoVector.segment(_alpha + a * _no_b, _no_b) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fia.segment(_alpha + a * _no_b, _no_b) * _trafoVector.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) + eigenvalue) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _trafoVector.segment(a * _no_a, _no_a) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _trafoVector.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) + eigenvalue) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getLeftAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getFockAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  // Fai is used to store -Fia (Right exc transformed Fock matrix)
  Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose() +
                        _gsLagrange.segment(i * _nv, _nv) * _Fai.segment(j * _nv, _nv).transpose() +
                        _Fai.segment(i * _nv, _nv) * _gsLagrange.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) - eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getFockAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getFockAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  // Fai is used to store -Fia (Right exc transformed Fock matrix)
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fai.segment(j * _nv_a, _nv_a).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) - eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Xia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _gsLagrange.segment(_alpha + i * _nv_b, _nv_b) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fai.segment(_alpha + i * _nv_b, _nv_b) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) - eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) - eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getFockAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getGLagrangeAmplitudes(unsigned i, unsigned j, int) {
  Eigen::MatrixXd amp = _Jia->middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Zia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Zia.middleRows(j * _nv, _nv).transpose() +
                        _gsLagrange.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _gsLagrange.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getGLagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getGLagrangeAmplitudes(unsigned i, unsigned j, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fia.segment(j * _nv_a, _nv_a).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j)) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Zia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _gsLagrange.segment(_alpha + i * _nv_b, _nv_b) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fia.segment(_alpha + i * _nv_b, _nv_b) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j)) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j)) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getGLagrangeAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getGLagrangeAmplitudesV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Jia->middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Zia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Zia.middleRows(b * _no, _no).transpose() +
                        _gsLagrange.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _gsLagrange.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getGLagrangeAmplitudesV() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getGLagrangeAmplitudesV(unsigned a, unsigned b, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Jia->alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _gsLagrange.segment(a * _no_a, _no_a) * _Fia.segment(b * _no_a, _no_a).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _gsLagrange.segment(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b)) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Zia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _gsLagrange.segment(_alpha + a * _no_b, _no_b) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fia.segment(_alpha + a * _no_b, _no_b) * _gsLagrange.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b)) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Jia->alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _gsLagrange.segment(a * _no_a, _no_a) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _gsLagrange.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b)) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getGLagrangeAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getELagrangeAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  // trafoVector is the excited state Lagrange multiplier in this case
  Eigen::MatrixXd amp = _Zia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Zia.middleRows(j * _nv, _nv).transpose() +
                        _trafoVector.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _trafoVector.segment(j * _nv, _nv).transpose() +
                        // Fock part of these amplitudes
                        _Xia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose() +
                        _gsLagrange.segment(i * _nv, _nv) * _Fai.segment(j * _nv, _nv).transpose() +
                        _Fai.segment(i * _nv, _nv) * _gsLagrange.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) - eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getELagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getELagrangeAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _trafoVector.segment(i * _nv_a, _nv_a) * _Fia.segment(j * _nv_a, _nv_a).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _trafoVector.segment(j * _nv_a, _nv_a).transpose() +
                          _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fai.segment(j * _nv_a, _nv_a).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) - eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _trafoVector.segment(_alpha + i * _nv_b, _nv_b) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fia.segment(_alpha + i * _nv_b, _nv_b) * _trafoVector.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Xia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _gsLagrange.segment(_alpha + i * _nv_b, _nv_b) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fai.segment(_alpha + i * _nv_b, _nv_b) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) - eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _trafoVector.segment(i * _nv_a, _nv_a) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _trafoVector.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _gsLagrange.segment(i * _nv_a, _nv_a) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) - eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getELagrangeAmplitudes() unrestricted */

template<>
Eigen::MatrixXd XWFController<RESTRICTED>::getELagrangeAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  // trafoVector is the excited state Lagrange multiplier in this case
  Eigen::MatrixXd amp = _Zia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Zia.middleRows(b * _no, _no).transpose() +
                        _trafoVector.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _trafoVector.segment(b * _no, _no).transpose() +
                        // Fock part of these amplitudes
                        _Xia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose() +
                        _gsLagrange.segment(a * _no, _no) * _Fai.segment(b * _no, _no).transpose() +
                        _Fai.segment(a * _no, _no) * _gsLagrange.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) - eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getELagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd XWFController<UNRESTRICTED>::getELagrangeAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _trafoVector.segment(a * _no_a, _no_a) * _Fia.segment(b * _no_a, _no_a).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _trafoVector.segment(b * _no_a, _no_a).transpose() +
                          _Xia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Xia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _gsLagrange.segment(a * _no_a, _no_a) * _Fai.segment(b * _no_a, _no_a).transpose() +
                          _Fai.segment(a * _no_a, _no_a) * _gsLagrange.segment(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) - eigenvalue) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _trafoVector.segment(_alpha + a * _no_b, _no_b) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fia.segment(_alpha + a * _no_b, _no_b) * _trafoVector.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Xia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _gsLagrange.segment(_alpha + a * _no_b, _no_b) * _Fai.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fai.segment(_alpha + a * _no_b, _no_b) * _gsLagrange.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) - eigenvalue) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _trafoVector.segment(a * _no_a, _no_a) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _trafoVector.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Xia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _gsLagrange.segment(a * _no_a, _no_a) * _Fai.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fai.segment(a * _no_a, _no_a) * _gsLagrange.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) - eigenvalue) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getELagrangeAmplitudesV() unrestricted */

template<Options::SCF_MODES SCFMode>
void XWFController<SCFMode>::reorderTensor(Eigen::Ref<Eigen::MatrixXd> Yia, unsigned n, unsigned m, unsigned l) {
  Timings::takeTime("Exc. State WF -   Tensor Reorder");
  for (size_t k = 0; k < l; ++k) {
    Eigen::Map<Eigen::MatrixXd> yian(Yia.data() + k * n * m, n, m);
    Eigen::Map<Eigen::MatrixXd> yiat(Yia.data() + k * n * m, m, n);
    yiat = yian.transpose().eval();
  }
  Timings::timeTaken("Exc. State WF -   Tensor Reorder");
} /* this->reorderTensor() */

template<Options::SCF_MODES SCFMode>
void XWFController<SCFMode>::performRITransformation(Eigen::MatrixXd& Yia, unsigned no, bool transposeRITrafo) {
  Timings::takeTime("Exc. State WF -     RI Transform");
  unsigned blockSize = Yia.rows() / no;
  bool naf = _nx < _nxb;

  if (transposeRITrafo && naf) {
    Yia.conservativeResize(Yia.rows(), _nxb);
  }

  for (size_t i = 0; i < no; ++i) {
    if (transposeRITrafo) {
      Yia.block(i * blockSize, 0, blockSize, _nxb) =
          (Yia.block(i * blockSize, 0, blockSize, _nx) * _riTrafo.transpose()).eval();
    }
    else {
      Yia.block(i * blockSize, 0, blockSize, _nx) = (Yia.block(i * blockSize, 0, blockSize, _nxb) * _riTrafo).eval();
    }
  }

  if (!transposeRITrafo && naf) {
    Yia.conservativeResize(Yia.rows(), _nx);
  }
  Timings::timeTaken("Exc. State WF -     RI Transform");
} /* this->performRITransformation() restricted/unrestricted */

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd XWFController<SCFMode>::getJ2GContribution(double* YiaPtr, double* JijPtr,
                                                           Eigen::Ref<Eigen::VectorXd> guessVector, bool fromLeft, int iSpin) {
  Timings::takeTime("Exc. State WF -     J2G Contrib.");
  // declarations
  unsigned _nv = (iSpin > 0) ? _nv_a : _nv_b;
  unsigned _no = (iSpin > 0) ? _no_a : _no_b;

  Eigen::MatrixXd coeffA;
  Eigen::MatrixXd coeffB;
  int fSpin = 1;
  for_spin(_P, _H) {
    if (iSpin * fSpin > 0) {
      coeffA = fromLeft ? _H_spin.middleCols(_no, _nv) : _P_spin.middleCols(_no, _nv);
      coeffB = fromLeft ? _P_spin.middleCols(_no, _nv) : _H_spin.middleCols(_no, _nv);
    }
    fSpin *= -1;
  };

  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nv * _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);

  unsigned nThreads = omp_get_max_threads();
  Eigen::MatrixXd Gia = Eigen::MatrixXd::Zero(_nb * _no, _nxb);
  Eigen::MatrixXd Eta = Eigen::MatrixXd::Zero(_nb * _no, nThreads);

  // contractions
  if (YiaPtr) {
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
      Eigen::Map<Eigen::MatrixXd> yia(YiaPtr + Q * _nv * _no, _nv, _no);
      gia.noalias() += coeffB * yia;
    }
  }
  if (JijPtr) {
    if (fromLeft) {
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
        Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
        gia.noalias() -= coeffB * (guess * jij.transpose());
      }
    }
    else {
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
        Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
        gia.noalias() -= coeffB * guess * jij;
      }
    }
  }

  Timings::timeTaken("Exc. State WF -     J2G Contrib.");
  this->performRITransformation(Gia, _no, true);
  Timings::takeTime("Exc. State WF -     J2G Contrib.");

  // loop over the rest
  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, size_t P, unsigned iThread) {
    Eigen::Map<Eigen::MatrixXd> eta(Eta.col(iThread).data(), _nb, _no);
    Eigen::Map<Eigen::MatrixXd> gia(Gia.col(P).data(), _nb, _no);
    eta.noalias() += AO * gia;
  };

  _integrals->loop(distribute);

  // sum over threads
  for (size_t iThread = 0; iThread < nThreads; ++iThread) {
    Eigen::Map<Eigen::MatrixXd> eta(Eta.col(iThread).data(), _nb, _no);
    sigma.noalias() += coeffA.transpose() * eta;
  }

  Timings::timeTaken("Exc. State WF -     J2G Contrib.");
  return sigmaVector;
} /* this->getJ2GContribution() restricted/unrestricted */

template<Options::SCF_MODES SCFMode>
void XWFController<SCFMode>::performTransformation(Eigen::MatrixXd& Xia, Eigen::Ref<Eigen::VectorXd> trafoVector,
                                                   bool fromLeft, int iSpin, bool singlesTransformed) {
  Timings::takeTime("Exc. State WF -    Int-Transform");

  // declarations
  unsigned _nv = (iSpin > 0) ? _nv_a : _nv_b;
  unsigned _no = (iSpin > 0) ? _no_a : _no_b;

  double* JijPtr;
  double* coeffAPtr;
  double* coeffBPtr;

  int fSpin = 1;
  auto& Jij = *_Jij;
  auto& Bij = *_Bij;
  for_spin(_C, _P, _H, Jij, Bij) {
    if (iSpin * fSpin > 0) {
      if (singlesTransformed) {
        coeffAPtr = fromLeft ? _H_spin.data() : _P_spin.data();
        coeffBPtr = fromLeft ? _P_spin.data() : _H_spin.data();
        JijPtr = Bij_spin.data();
      }
      else {
        coeffAPtr = _C_spin.data();
        coeffBPtr = _C_spin.data();
        JijPtr = Jij_spin.data();
      }
    }
    fSpin *= -1;
  };

  Eigen::Map<Eigen::MatrixXd> coeffA(coeffAPtr + _nb * _no, _nb, _nv);
  Eigen::Map<Eigen::MatrixXd> coeffB(coeffBPtr + _nb * _no, _nb, _nv);

  Eigen::Map<Eigen::MatrixXd> trafo(trafoVector.data(), _nv, _no);
  Eigen::MatrixXd eta = coeffB * trafo;

  Xia.resize(_nb * _no, _nxb);

  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, unsigned long P, unsigned) {
    Eigen::Map<Eigen::MatrixXd> xia(Xia.col(P).data(), _nb, _no);
    xia.noalias() = AO * eta;
  };

  _integrals->loop(distribute);
  this->performRITransformation(Xia, _no);

  // Jij contribution
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> xia_b(Xia.col(Q).data(), _nb, _no);
    Eigen::Map<Eigen::MatrixXd> xia_v(Xia.col(Q).data(), _nv, _no);
    xia_v = coeffA.rightCols(_nv).transpose() * xia_b;

    if (fromLeft) {
      xia_v.noalias() -= trafo * jij.transpose();
    }
    else {
      xia_v.noalias() -= trafo * jij;
    }
  }

  Xia.conservativeResize(_nv * _no, _nx);
  Timings::timeTaken("Exc. State WF -    Int-Transform");
} /* this->performTransformation() restricted/unrestricted */

template<>
void XWFController<RESTRICTED>::transformIntegrals() {
  unsigned mo = _no + _nv;
  Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(mo, mo);
  t1.bottomLeftCorner(_nv, _no) = Eigen::Map<Eigen::MatrixXd>(_singles.data(), _nv, _no);
  _P = _C;
  _H = _C;
  _P.leftCols(mo) = _C.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) - t1.transpose());
  _H.leftCols(mo) = _C.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) + t1);

  // here I'm misusing the transformation function since basically everything it does is:
  // _Bia->col(Q) = jab(Q) * sgl - sgl * jij(Q), which just happens to be part of Bia
  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);
  this->performTransformation(*_Bia, _singles, false, 1, false);

  (*_Bij) = (*_Jij);
  (*_Bia) += (*_Jia);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> bia(_Bia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    bij.noalias() += jia.transpose() * sgl;
    bia.noalias() -= sgl * jia.transpose() * sgl;
  }
} /* this->transformIntegrals() */

template<>
void XWFController<UNRESTRICTED>::transformIntegrals() {
  unsigned iSpin = 0;
  for_spin(_P, _H, _C, _nv, _no) {
    unsigned mo = _no_spin + _nv_spin;
    Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(mo, mo);
    t1.bottomLeftCorner(_nv_spin, _no_spin) =
        Eigen::Map<Eigen::MatrixXd>(_singles.data() + iSpin * _alpha, _nv_spin, _no_spin);
    _P_spin = _C_spin;
    _H_spin = _C_spin;
    _P_spin.leftCols(mo) = _C_spin.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) - t1.transpose());
    _H_spin.leftCols(mo) = _C_spin.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) + t1);
    ++iSpin;
  };

  // here I'm misusing the transformation function since basically everything it does is:
  // _Bia->col(Q) = jab(Q) * sgl - sgl * jij(Q), which just happens to be part of Bia
  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);
  this->performTransformation(_Bia->alpha, _singles.head(_alpha), false, 1, false);
  this->performTransformation(_Bia->beta, _singles.tail(_beta), false, -1, false);
  _Bij->alpha = _Jij->alpha;
  _Bia->alpha += _Jia->alpha;
  _Bij->beta = _Jij->beta;
  _Bia->beta += _Jia->beta;
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> bia_a(_Bia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    bij_a.noalias() += jia_a.transpose() * sgl_a;
    bia_a.noalias() -= sgl_a * jia_a.transpose() * sgl_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> bia_b(_Bia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    bij_b.noalias() += jia_b.transpose() * sgl_b;
    bia_b.noalias() -= sgl_b * jia_b.transpose() * sgl_b;
  }
} /* this->transformIntegrals() */

template<>
void XWFController<RESTRICTED>::calculateGroundstate() {
  // matrices to store vectors for diis
  Eigen::MatrixXd sglStorage, resStorage;

  // timings
  std::chrono::steady_clock::time_point itStart = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point itEnd = std::chrono::steady_clock::now();

  printf("    Convergence threshold             : %-5.1e\n", _settings.conv);
  printf("    Maximum number of iterations      : %-5i\n\n", _settings.maxCycles);
  printf("    The zeroth iteration corresponds to MP2.\n\n");

  _iter = 0;
  printf("   it.   time (min)      gs energy      corr energy      res norm      t1 norm  \n");
  printf("  ------------------------------------------------------------------------------\n");
  do {
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // update singles with residuals
      for (size_t i = 0; i < _no; ++i) {
        for (size_t a = 0; a < _nv; ++a) {
          _singles(i * _nv + a) += _residual(i * _nv + a) / (_e(i) - _e(_no + a));
        }
      }
      if (_iter > 0) {
        if (_iter <= _settings.diisStore) {
          resStorage.conservativeResize(_no * _nv, resStorage.cols() + 1);
          sglStorage.conservativeResize(_no * _nv, sglStorage.cols() + 1);
        }
        resStorage.col((_iter - 1) % _settings.diisStore) = _residual;
        sglStorage.col((_iter - 1) % _settings.diisStore) = _singles;

        Eigen::VectorXd rhs(Eigen::VectorXd::Zero(_iter + 1));
        Eigen::MatrixXd lhs(Eigen::MatrixXd::Zero(_iter + 1, _iter + 1));

        for (size_t iStored = 0; iStored < _iter; ++iStored) {
          lhs(iStored, _iter) = -1.0;
          lhs(_iter, iStored) = -1.0;
        }
        lhs.topLeftCorner(_iter, _iter) = resStorage.transpose() * resStorage;
        double norm = lhs.topLeftCorner(_iter, _iter).cwiseAbs().maxCoeff();
        lhs.topLeftCorner(_iter, _iter) *= 1.0 / norm;
        rhs(_iter) = -1.0;
        _singles = sglStorage * lhs.householderQr().solve(rhs).head(_iter);
      }
      this->transformIntegrals();
      this->calculateResidual();
    }
    else {
      // Yia only needs to be calculated once
      _Yia.setZero();
      if (_settings.ltconv == 0) {
        for (size_t i = 0; i < _no; ++i) {
          for (size_t j = i; j < _no; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j);
            _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Jia->middleRows(j * _nv, _nv);
            if (i != j) {
              _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Jia->middleRows(i * _nv, _nv);
            }
          }
        }
      }
      else {
        for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
          Eigen::VectorXd param(_nv * _no);
          for (size_t i = 0, ia = 0; i < _no; ++i) {
            for (size_t a = _no; a < _no + _nv; ++a, ++ia) {
              param(ia) = std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
            }
          }
          Eigen::MatrixXd QP = (*_Jia).transpose() * param.asDiagonal() * (*_Jia);
          _Yia.noalias() -= _oss * _weights(iRoot) * param.asDiagonal() * (*_Jia) * QP;
        }
      }
      _P = _C;
      _H = _C;
    }

    // calculate correlation energy
    _corrEnergy = _Yia.cwiseProduct(*_Jia).sum() + _Fia.dot(_singles);
    if (_iter == 0) {
      _mp2Energy = _corrEnergy;
    }

    // iteration time
    itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();

    // print info
    printf("%5i %11.3f %18.10f %15.10f %14.10f %10.6f\n", _iter++, duration / 60.0, _hfEnergy + _corrEnergy,
           _corrEnergy, _residual.norm(), _singles.norm());

    itStart = std::chrono::steady_clock::now();

    if (_iter > _settings.maxCycles) {
      printf("\n  Could not converge t1 amplitudes.\n");
      break;
    }
  } while (_residual.norm() > _settings.conv && _xwfModel == Options::LR_METHOD::CC2);
  printf("\n    RHF energy (a.u.): %24.10f\n", _hfEnergy);
  printf("    Final MP2 energy (a.u.): %18.10f\n", _mp2Energy);
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printf("    Final CC2 energy (a.u.): %18.10f\n", _corrEnergy);
  }
  printf("\n");
} /* this->calculateGroundstate() restricted */

template<>
void XWFController<UNRESTRICTED>::calculateGroundstate() {
  // matrices to store vectors for diis
  Eigen::MatrixXd sglStorage, resStorage;

  // timings
  std::chrono::steady_clock::time_point itStart = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point itEnd = std::chrono::steady_clock::now();

  printf("    Convergence threshold             : %-5.1e\n", _settings.conv);
  printf("    Maximum number of iterations      : %-5i\n\n", _settings.maxCycles);
  printf("    The zeroth iteration corresponds to MP2.\n\n");

  _iter = 0;
  printf("   it.   time (min)      gs energy      corr energy      res norm      t1 norm  \n");
  printf("  ------------------------------------------------------------------------------\n");
  do {
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // update singles with residuals
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t a = 0; a < _nv_a; ++a) {
          _singles(i * _nv_a + a) += _residual(i * _nv_a + a) / (_e.alpha(i) - _e.alpha(_no_a + a));
        }
      }
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t a = 0; a < _nv_b; ++a) {
          _singles(_alpha + i * _nv_b + a) += _residual(_alpha + i * _nv_b + a) / (_e.beta(i) - _e.beta(_no_b + a));
        }
      }
      if (_iter > 0) {
        if (_iter <= _settings.diisStore) {
          resStorage.conservativeResize(_nDim, resStorage.cols() + 1);
          sglStorage.conservativeResize(_nDim, sglStorage.cols() + 1);
        }
        resStorage.col((_iter - 1) % _settings.diisStore) = _residual;
        sglStorage.col((_iter - 1) % _settings.diisStore) = _singles;

        Eigen::VectorXd rhs(Eigen::VectorXd::Zero(_iter + 1));
        Eigen::MatrixXd lhs(Eigen::MatrixXd::Zero(_iter + 1, _iter + 1));

        for (size_t iStored = 0; iStored < _iter; ++iStored) {
          lhs(iStored, _iter) = -1.0;
          lhs(_iter, iStored) = -1.0;
        }
        lhs.topLeftCorner(_iter, _iter) = resStorage.transpose() * resStorage;
        double norm = lhs.topLeftCorner(_iter, _iter).cwiseAbs().maxCoeff();
        lhs.topLeftCorner(_iter, _iter) *= 1.0 / norm;
        rhs(_iter) = -1.0;
        _singles = sglStorage * lhs.householderQr().solve(rhs).head(_iter);
      }
      this->transformIntegrals();
      this->calculateResidual();
    }
    else {
      if (_settings.ltconv == 0) {
        // alpha
        for (size_t i = 0; i < _no_a; ++i) {
          for (size_t j = i; j < _no_a; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
            _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
            if (i != j) {
              _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
            }
          }
        }

        // beta
        for (size_t i = 0; i < _no_b; ++i) {
          for (size_t j = i; j < _no_b; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
            _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
            if (i != j) {
              _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
            }
          }
        }

        // mixed
        for (size_t i = 0; i < _no_a; ++i) {
          for (size_t j = 0; j < _no_b; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
            _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
            _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
          }
        }
      }
      else {
        for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
          Eigen::VectorXd param_a(_nv_a * _no_a);
          for (size_t i = 0, ia = 0; i < _no_a; ++i) {
            for (size_t a = _no_a; a < _no_a + _nv_a; ++a, ++ia) {
              param_a(ia) = std::exp(-(_e.alpha(a) - _e.alpha(i)) * _roots(iRoot));
            }
          }
          Eigen::VectorXd param_b(_nv_b * _no_b);
          for (size_t j = 0, jb = 0; j < _no_b; ++j) {
            for (size_t b = _no_b; b < _no_b + _nv_b; ++b, ++jb) {
              param_b(jb) = std::exp(-(_e.beta(b) - _e.beta(j)) * _roots(iRoot));
            }
          }
          Eigen::MatrixXd QP_a = _Jia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha;
          Eigen::MatrixXd QP_b = _Jia->beta.transpose() * param_b.asDiagonal() * _Jia->beta;

          _Yia.alpha.noalias() -= _oss * _weights(iRoot) * param_a.asDiagonal() * _Jia->alpha * QP_b;
          _Yia.beta.noalias() -= _oss * _weights(iRoot) * param_b.asDiagonal() * _Jia->beta * QP_a;
        }
      }
      for_spin(_C, _P, _H) {
        _P_spin = _C_spin;
        _H_spin = _C_spin;
      };
    }

    // calculate correlation energy
    _corrEnergy = 0.5 * _Fia.dot(_singles);
    auto& Jia = *_Jia;
    for_spin(_Yia, Jia) {
      _corrEnergy += 0.5 * _Yia_spin.cwiseProduct(Jia_spin).sum();
    };
    if (_iter == 0) {
      _mp2Energy = _corrEnergy;
    }

    // iteration time
    itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();

    // print info
    printf("%5i %11.3f %18.10f %15.10f %14.10f %10.6f\n", _iter, duration / 60.0, _hfEnergy + _corrEnergy, _corrEnergy,
           _residual.norm(), _singles.norm());

    itStart = std::chrono::steady_clock::now();

    if (_iter > _settings.maxCycles) {
      printf("\n  Could not converge t1 amplitudes.\n");
      break;
    }
    ++_iter;
  } while (_residual.norm() > _settings.conv && _xwfModel == Options::LR_METHOD::CC2);
  printf("\n    UHF energy (a.u.): %24.10f\n", _hfEnergy);
  printf("    Final MP2 energy (a.u.): %18.10f\n", _mp2Energy);
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printf("    Final CC2 energy (a.u.): %18.10f\n", _corrEnergy);
  }
  printf("\n");
} /* this->calculateGroundstate() unrestricted */

template<>
void XWFController<RESTRICTED>::initialize() {
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

  _Yia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);

  _Fai = Eigen::VectorXd::Zero(_no * _nv);
  _Fia = Eigen::VectorXd::Zero(_no * _nv);
  _singles = Eigen::VectorXd::Zero(_no * _nv);
  _residual = Eigen::VectorXd::Zero(_no * _nv);
  _gsLagrange = Eigen::VectorXd::Zero(_no * _nv);

  _trafoVector = Eigen::VectorXd::Zero(_no * _nv);

  _eVirt = Eigen::MatrixXd::Zero(_nv, _nv);
  _eVirt.colwise() += _e.segment(_no, _nv);
  _eVirt.rowwise() += _e.segment(_no, _nv).transpose();

  _eOcc = Eigen::MatrixXd::Zero(_no, _no);
  _eOcc.colwise() += _e.segment(0, _no);
  _eOcc.rowwise() += _e.segment(0, _no).transpose();

  Timings::takeTime("Exc. State WF -     Ground State");
  this->calculateGroundstate();
  Timings::timeTaken("Exc. State WF -     Ground State");

  if (_settings.nEigen > 0) {
    Timings::takeTime("Exc. State WF -   E-Intermediate");
    this->calculateE();
    Timings::timeTaken("Exc. State WF -   E-Intermediate");
  }
} /* this->initialize() restricted */

template<>
void XWFController<UNRESTRICTED>::initialize() {
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

  _Yia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  for_spin(_Yia, _nv, _no) {
    _Yia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
  };

  _Fai = Eigen::VectorXd::Zero(_nDim);
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _singles = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);
  _gsLagrange = Eigen::VectorXd::Zero(_nDim);

  _trafoVector = Eigen::VectorXd::Zero(_nDim);

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

  Timings::takeTime("Exc. State WF -     Ground State");
  this->calculateGroundstate();
  Timings::timeTaken("Exc. State WF -     Ground State");

  if (_settings.nEigen > 0) {
    Timings::takeTime("Exc. State WF -   E-Intermediate");
    this->calculateE();
    Timings::timeTaken("Exc. State WF -   E-Intermediate");
  }
} /* this->initialize() unrestricted */

template<>
void XWFController<RESTRICTED>::normalizeEigenvectors(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues) {
  Timings::takeTime("Exc. State WF -    Normalization");

  if (eigenvectors.size() == 1) {
    printf("  Normalizing eigenvectors such that <~R|R> = 1.\n\n");
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      double crr = eigenvector.dot(eigenvector);

      this->performTransformation(_Xia, eigenvector, false);
      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
          crr += (i == j ? 0.5 : 1.0) * (_soss * rij - _sss * rij.transpose()).cwiseProduct(rij).sum();
        }
      }
      eigenvector *= 1 / std::sqrt(crr);
    } /* iEigen */
    eigenvectors.push_back(eigenvectors[0]);
  }
  else {
    printf("  Normalizing eigenvectors such that <L|R> = <R|R> = 1.\n\n");
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
      Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      double clr = leftEigenvector.dot(rightEigenvector);
      double crr = rightEigenvector.dot(rightEigenvector);
      this->performTransformation(_Xia, rightEigenvector, false);
      this->performTransformation(_Zia, leftEigenvector, true);
      _trafoVector = leftEigenvector;

      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }
      rightEigenvector *= 1 / std::sqrt(crr);
      leftEigenvector *= std::sqrt(crr) / clr;
    } /* iEigen */
  }
  Timings::timeTaken("Exc. State WF -    Normalization");
} /* this->normalizeEigenvectors() restricted */

template<>
void XWFController<UNRESTRICTED>::normalizeEigenvectors(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues) {
  Timings::takeTime("Exc. State WF -    Normalization");

  if (eigenvectors.size() == 1) {
    printf("  Normalizing eigenvectors such that <~R|R> = 1.\n\n");
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      // get eigenpair
      Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      // calculate scaling factor
      double crr = eigenvector.dot(eigenvector);
      this->performTransformation(_Xia.alpha, eigenvector.head(_alpha), false, 1);
      this->performTransformation(_Xia.beta, eigenvector.tail(_beta), false, -1);

      // alpha
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
          crr += (i == j ? 0.5 : 1.0) * (_sss * rij - _sss * rij.transpose()).cwiseProduct(rij).sum();
        }
      }

      // beta
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
          crr += (i == j ? 0.5 : 1.0) * (_sss * rij - _sss * rij.transpose()).cwiseProduct(rij).sum();
        }
      }

      // mixed (double counting for alpha -> beta and beta -> alpha)
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
          crr += (0.5 + 0.5) * (_oss * rij).cwiseProduct(rij).sum();
        }
      }

      // scale eigenvector
      eigenvector *= 1 / std::sqrt(crr);
    } /* iEigen */
    eigenvectors.push_back(eigenvectors[0]);
  }
  else {
    printf("  Normalizing eigenvectors such that <L|R> = <R|R> = 1.\n\n");
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
      Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      // calculate scaling factor
      double clr = leftEigenvector.dot(rightEigenvector);
      double crr = rightEigenvector.dot(rightEigenvector);
      this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
      this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
      this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
      this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);
      _trafoVector = leftEigenvector;

      // alpha
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }

      // beta
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }

      // mixed (double counting for alpha -> beta and beta -> alpha)
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);
          clr += (0.5 + 0.5) * lij.cwiseProduct(rij).sum();
          crr += (0.5 + 0.5) * rij.cwiseProduct(rij).sum();
        }
      }

      rightEigenvector *= 1 / std::sqrt(crr);
      leftEigenvector *= std::sqrt(crr) / clr;
    } /* iEigen */
  }
  Timings::timeTaken("Exc. State WF -    Normalization");
} /* this->normalizeEigenvectors() unrestricted */

template class XWFController<Options::SCF_MODES::RESTRICTED>;
template class XWFController<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity