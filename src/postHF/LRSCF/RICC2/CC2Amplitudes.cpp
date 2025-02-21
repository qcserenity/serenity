/**
 * @file CC2Amplitudes.cpp
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
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getAmplitudes(unsigned i, unsigned j, int) {
  // Be careful: in the triplet case, the MP2 calculation and the calculation of the Yia
  // intermediate do not use (_soss = _sss - _oss) so I don't have an if-else block here.
  // (_soss = _sss - _oss) is being set after calculation of the E-intermediate!
  Eigen::MatrixXd amp = _Bia->middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getAmplitudes(unsigned i, unsigned j, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getAmplitudesA(unsigned i, unsigned j, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getAmplitudesA() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getAmplitudesA(unsigned i, unsigned j, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getAmplitudesV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getAmplitudesV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getAmplitudesV(unsigned a, unsigned b, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getAmplitudesAV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Bia->middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return amp.cwiseQuotient(denom);
} /* this->getAmplitudesAV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getAmplitudesAV(unsigned a, unsigned b, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getRightAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  if (_settings.triplet) {
    Eigen::MatrixXd amp1 = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
    Eigen::MatrixXd amp2 = _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
    return ((_sss + _oss) * amp1 + (_sss - _oss) * amp2 - _sss * (amp1 + amp2).transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose() +
                          _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
    return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
} /* this->getRightAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getRightAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getRightAmplitudesA(unsigned i, unsigned j, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose() +
                        _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getRightAmplitudesA() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getRightAmplitudesA(unsigned i, unsigned j, double eigenvalue, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getRightAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose() +
                        _Bia->middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getRightAmplitudesV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getRightAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getRightAmplitudesAV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(a * _no, _no) * _Bia->middleRows(b * _no, _no).transpose() +
                        _Bia->middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return amp.cwiseQuotient(denom);
} /* this->getRightAmplitudesAV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getRightAmplitudesAV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getLeftAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Zia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Zia.middleRows(j * _nv, _nv).transpose() +
                        _leftSingles.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _leftSingles.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getLeftAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getLeftAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _leftSingles.segment(i * _nv_a, _nv_a) * _Fia.segment(j * _nv_a, _nv_a).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _leftSingles.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) + eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _leftSingles.segment(_alpha + i * _nv_b, _nv_b) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fia.segment(_alpha + i * _nv_b, _nv_b) * _leftSingles.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) + eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Zia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _leftSingles.segment(i * _nv_a, _nv_a) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _leftSingles.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) + eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getLeftAmplitudes() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getLeftAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  Eigen::MatrixXd amp = _Zia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Zia.middleRows(b * _no, _no).transpose() +
                        _leftSingles.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _leftSingles.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getLeftAmplitudesV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getLeftAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _leftSingles.segment(a * _no_a, _no_a) * _Fia.segment(b * _no_a, _no_a).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _leftSingles.segment(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) + eigenvalue) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Zia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _leftSingles.segment(_alpha + a * _no_b, _no_b) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fia.segment(_alpha + a * _no_b, _no_b) * _leftSingles.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) + eigenvalue) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Zia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Zia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _leftSingles.segment(a * _no_a, _no_a) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _leftSingles.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) + eigenvalue) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getLeftAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getFockAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  // Fai is used to store -Fia (Right exc transformed Fock matrix)
  Eigen::MatrixXd amp = _Wia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Wia.middleRows(j * _nv, _nv).transpose() +
                        _fockSingles.segment(i * _nv, _nv) * _Fai.segment(j * _nv, _nv).transpose() +
                        _Fai.segment(i * _nv, _nv) * _fockSingles.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) - eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getFockAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getFockAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  // Fai is used to store -Fia (Right exc transformed Fock matrix)
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Wia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Wia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _fockSingles.segment(i * _nv_a, _nv_a) * _Fai.segment(j * _nv_a, _nv_a).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _fockSingles.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) - eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Wia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Wia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _fockSingles.segment(_alpha + i * _nv_b, _nv_b) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fai.segment(_alpha + i * _nv_b, _nv_b) * _fockSingles.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) - eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Wia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Wia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _fockSingles.segment(i * _nv_a, _nv_a) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _fockSingles.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) - eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getFockAmplitudes() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getGLagrangeAmplitudes(unsigned i, unsigned j, int) {
  Eigen::MatrixXd amp = _Jia->middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Zia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Zia.middleRows(j * _nv, _nv).transpose() +
                        _gsLagrange.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _gsLagrange.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getGLagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getGLagrangeAmplitudes(unsigned i, unsigned j, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getGLagrangeAmplitudesV(unsigned a, unsigned b, int) {
  Eigen::MatrixXd amp = _Jia->middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Zia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Zia.middleRows(b * _no, _no).transpose() +
                        _gsLagrange.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _gsLagrange.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b)) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getGLagrangeAmplitudesV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getGLagrangeAmplitudesV(unsigned a, unsigned b, int iSpin) {
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
Eigen::MatrixXd CC2Controller<RESTRICTED>::getELagrangeAmplitudes(unsigned i, unsigned j, double eigenvalue, int) {
  // trafoVector is the excited-state Lagrange multiplier in this case
  Eigen::MatrixXd amp = _Yia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Yia.middleRows(j * _nv, _nv).transpose() +
                        _exlSingles1.segment(i * _nv, _nv) * _Fia.segment(j * _nv, _nv).transpose() +
                        _Fia.segment(i * _nv, _nv) * _exlSingles1.segment(j * _nv, _nv).transpose() +
                        // Fock part of these amplitudes
                        _Wia.middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose() +
                        _Jia->middleRows(i * _nv, _nv) * _Wia.middleRows(j * _nv, _nv).transpose() +
                        _exlSingles2.segment(i * _nv, _nv) * _Fai.segment(j * _nv, _nv).transpose() +
                        _Fai.segment(i * _nv, _nv) * _exlSingles2.segment(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) - eigenvalue) - _eVirt;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getELagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getELagrangeAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Yia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Yia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _exlSingles1.segment(i * _nv_a, _nv_a) * _Fia.segment(j * _nv_a, _nv_a).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _exlSingles1.segment(j * _nv_a, _nv_a).transpose() +
                          _Wia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Wia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _exlSingles2.segment(i * _nv_a, _nv_a) * _Fai.segment(j * _nv_a, _nv_a).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _exlSingles2.segment(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) - eigenvalue) - _eVirt.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Yia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Yia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _exlSingles1.segment(_alpha + i * _nv_b, _nv_b) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fia.segment(_alpha + i * _nv_b, _nv_b) * _exlSingles1.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Wia.beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Wia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
        _exlSingles2.segment(_alpha + i * _nv_b, _nv_b) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
        _Fai.segment(_alpha + i * _nv_b, _nv_b) * _exlSingles2.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) - eigenvalue) - _eVirt.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Yia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Yia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _exlSingles1.segment(i * _nv_a, _nv_a) * _Fia.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fia.segment(i * _nv_a, _nv_a) * _exlSingles1.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Wia.alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Wia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _exlSingles2.segment(i * _nv_a, _nv_a) * _Fai.segment(_alpha + j * _nv_b, _nv_b).transpose() +
                          _Fai.segment(i * _nv_a, _nv_a) * _exlSingles2.segment(_alpha + j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) - eigenvalue) - _eVirtab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getELagrangeAmplitudes() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getELagrangeAmplitudesV(unsigned a, unsigned b, double eigenvalue, int) {
  // trafoVector is the excited-state Lagrange multiplier in this case
  Eigen::MatrixXd amp = _Yia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Yia.middleRows(b * _no, _no).transpose() +
                        _exlSingles1.segment(a * _no, _no) * _Fia.segment(b * _no, _no).transpose() +
                        _Fia.segment(a * _no, _no) * _exlSingles1.segment(b * _no, _no).transpose() +
                        // Fock part of these amplitudes
                        _Wia.middleRows(a * _no, _no) * _Jia->middleRows(b * _no, _no).transpose() +
                        _Jia->middleRows(a * _no, _no) * _Wia.middleRows(b * _no, _no).transpose() +
                        _exlSingles2.segment(a * _no, _no) * _Fai.segment(b * _no, _no).transpose() +
                        _Fai.segment(a * _no, _no) * _exlSingles2.segment(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) - eigenvalue) + _eOcc;
  return (_soss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
} /* this->getELagrangeAmplitudes() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getELagrangeAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Yia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Yia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _exlSingles1.segment(a * _no_a, _no_a) * _Fia.segment(b * _no_a, _no_a).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _exlSingles1.segment(b * _no_a, _no_a).transpose() +
                          _Wia.alpha.middleRows(a * _no_a, _no_a) * _Jia->alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Wia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _exlSingles2.segment(a * _no_a, _no_a) * _Fai.segment(b * _no_a, _no_a).transpose() +
                          _Fai.segment(a * _no_a, _no_a) * _exlSingles2.segment(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) - eigenvalue) + _eOcc.alpha;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp =
        _Yia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Yia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _exlSingles1.segment(_alpha + a * _no_b, _no_b) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fia.segment(_alpha + a * _no_b, _no_b) * _exlSingles1.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Wia.beta.middleRows(a * _no_b, _no_b) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
        _Jia->beta.middleRows(a * _no_b, _no_b) * _Wia.beta.middleRows(b * _no_b, _no_b).transpose() +
        _exlSingles2.segment(_alpha + a * _no_b, _no_b) * _Fai.segment(_alpha + b * _no_b, _no_b).transpose() +
        _Fai.segment(_alpha + a * _no_b, _no_b) * _exlSingles2.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) - eigenvalue) + _eOcc.beta;
    return (_sss * amp - _sss * amp.transpose()).cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Yia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Yia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _exlSingles1.segment(a * _no_a, _no_a) * _Fia.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fia.segment(a * _no_a, _no_a) * _exlSingles1.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Wia.alpha.middleRows(a * _no_a, _no_a) * _Jia->beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Jia->alpha.middleRows(a * _no_a, _no_a) * _Wia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _exlSingles2.segment(a * _no_a, _no_a) * _Fai.segment(_alpha + b * _no_b, _no_b).transpose() +
                          _Fai.segment(a * _no_a, _no_a) * _exlSingles2.segment(_alpha + b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) - eigenvalue) + _eOccab;
    return (_oss * amp).cwiseQuotient(denom);
  }
} /* this->getELagrangeAmplitudesV() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getXiDoubles(unsigned i, unsigned j, double frequency,
                                                        SPMatrix<RESTRICTED>& dipoles, int) {
  Eigen::MatrixXd ij = this->getAmplitudes(i, j);
  Eigen::MatrixXd amp = ij * dipoles.transpose() + dipoles * ij;
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + frequency) - _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getXiDoubles() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getXiDoubles(unsigned i, unsigned j, double frequency,
                                                          SPMatrix<UNRESTRICTED>& dipoles, int iSpin) {
  Eigen::MatrixXd ij = this->getAmplitudes(i, j, iSpin);
  if (iSpin > 0) {
    Eigen::MatrixXd amp = ij * dipoles.alpha.transpose() + dipoles.alpha * ij;
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j) + frequency) - _eVirt.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = ij * dipoles.beta.transpose() + dipoles.beta * ij;
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j) + frequency) - _eVirt.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = ij * dipoles.beta.transpose() + dipoles.alpha * ij;
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j) + frequency) - _eVirtab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getXiDoubles() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getXiDoublesV(unsigned a, unsigned b, double frequency,
                                                         SPMatrix<RESTRICTED>& dipoles, int) {
  Eigen::MatrixXd ab = this->getAmplitudesV(a, b);
  Eigen::MatrixXd amp = ab * dipoles + dipoles.transpose() * ab;
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, -_e(_no + a) - _e(_no + b) + frequency) + _eOcc;
  return -1.0 * amp.cwiseQuotient(denom);
} /* this->getXiDoublesV() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getXiDoublesV(unsigned a, unsigned b, double frequency,
                                                           SPMatrix<UNRESTRICTED>& dipoles, int iSpin) {
  Eigen::MatrixXd ab = this->getAmplitudesV(a, b, iSpin);
  if (iSpin > 0) {
    Eigen::MatrixXd amp = ab * dipoles.alpha + dipoles.alpha.transpose() * ab;
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_a, -_e.alpha(_no_a + a) - _e.alpha(_no_a + b) + frequency) + _eOcc.alpha;
    return -1.0 * amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = ab * dipoles.beta + dipoles.beta.transpose() * ab;
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_b, _no_b, -_e.beta(_no_b + a) - _e.beta(_no_b + b) + frequency) + _eOcc.beta;
    return -1.0 * amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = ab * dipoles.beta + dipoles.alpha.transpose() * ab;
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_no_a, _no_b, -_e.alpha(_no_a + a) - _e.beta(_no_b + b) + frequency) + _eOccab;
    return -1.0 * amp.cwiseQuotient(denom);
  }
} /* this->getXiDoublesV() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getLTDoubles(unsigned i, unsigned j, double frequency, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(i * _nv, _nv) * _Yia.middleRows(j * _nv, _nv).transpose() +
                        _Yia.middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, -_e(i) - _e(j) - frequency) + _eVirt;
  return amp.cwiseQuotient(denom);
} /* this->getLTDoubles() restricted */

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::getLTDoublesV(unsigned a, unsigned b, double frequency, int) {
  Eigen::MatrixXd amp = _Xia.middleRows(a * _no, _no) * _Yia.middleRows(b * _no, _no).transpose() +
                        _Yia.middleRows(a * _no, _no) * _Xia.middleRows(b * _no, _no).transpose();
  Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_no, _no, _e(_no + a) + _e(_no + b) - frequency) - _eOcc;
  return amp.cwiseQuotient(denom);
} /* this->getLTDoublesV() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getLTDoubles(unsigned i, unsigned j, double frequency, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Yia.alpha.middleRows(j * _nv_a, _nv_a).transpose() +
                          _Yia.alpha.middleRows(i * _nv_a, _nv_a) * _Xia.alpha.middleRows(j * _nv_a, _nv_a).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, -_e.alpha(i) - _e.alpha(j) - frequency) + _eVirt.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(i * _nv_b, _nv_b) * _Yia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Yia.beta.middleRows(i * _nv_b, _nv_b) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, -_e.beta(i) - _e.beta(j) - frequency) + _eVirt.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(i * _nv_a, _nv_a) * _Yia.beta.middleRows(j * _nv_b, _nv_b).transpose() +
                          _Yia.alpha.middleRows(i * _nv_a, _nv_a) * _Xia.beta.middleRows(j * _nv_b, _nv_b).transpose();
    Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, -_e.alpha(i) - _e.beta(j) - frequency) + _eVirtab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getLTDoubles() unrestricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::getLTDoublesV(unsigned a, unsigned b, double frequency, int iSpin) {
  if (iSpin > 0) {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Yia.alpha.middleRows(b * _no_a, _no_a).transpose() +
                          _Yia.alpha.middleRows(a * _no_a, _no_a) * _Xia.alpha.middleRows(b * _no_a, _no_a).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(_no_a + a) + _e.alpha(_no_a + b) - frequency) - _eOcc.alpha;
    return amp.cwiseQuotient(denom);
  }
  else if (iSpin < 0) {
    Eigen::MatrixXd amp = _Xia.beta.middleRows(a * _no_b, _no_b) * _Yia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Yia.beta.middleRows(a * _no_b, _no_b) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(_no_b + a) + _e.beta(_no_b + b) - frequency) - _eOcc.beta;
    return amp.cwiseQuotient(denom);
  }
  else {
    Eigen::MatrixXd amp = _Xia.alpha.middleRows(a * _no_a, _no_a) * _Yia.beta.middleRows(b * _no_b, _no_b).transpose() +
                          _Yia.alpha.middleRows(a * _no_a, _no_a) * _Xia.beta.middleRows(b * _no_b, _no_b).transpose();
    Eigen::MatrixXd denom =
        Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(_no_a + a) + _e.beta(_no_b + b) - frequency) - _eOccab;
    return amp.cwiseQuotient(denom);
  }
} /* this->getLTDoublesV() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity