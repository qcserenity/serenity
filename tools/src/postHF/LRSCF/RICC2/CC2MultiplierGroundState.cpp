/**
 * @file CC2MultiplierGroundState.cpp
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

template<Options::SCF_MODES SCFMode>
void CC2Controller<SCFMode>::setGroundStateLagrangeMultipler(Eigen::Ref<Eigen::VectorXd> gsLagrange) {
  _gsLagrange = gsLagrange;
}

template<>
Eigen::VectorXd CC2Controller<RESTRICTED>::calculateGroundStateLagrangeMultiplier() {
  _Zia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd iajb = _Jia->middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
      Eigen::MatrixXd amps = (_soss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.middleRows(i * _nv, _nv).noalias() += amps * _Bia->middleRows(j * _nv, _nv);
      if (i != j) {
        _Zia.middleRows(j * _nv, _nv).noalias() += amps.transpose() * _Bia->middleRows(i * _nv, _nv);
      }
    }
  }

  Eigen::MatrixXd rightHandSide = _Fia;
  Eigen::Map<Eigen::MatrixXd> rhs(rightHandSide.data(), _nv, _no);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
    rhs.noalias() -= zia * bij.transpose();
  }
  Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nv * _no);
  rightHandSide.col(0) += this->getJ2GContribution(_Zia.data(), nullptr, nullVector, true);
  rightHandSide *= -1;

  return rightHandSide;
} /* this->calculateGroundStateLagrangeMultiplier() restricted */

template<>
Eigen::VectorXd CC2Controller<UNRESTRICTED>::calculateGroundStateLagrangeMultiplier() {
  _Zia.alpha.setZero();
  _Zia.beta.setZero();
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j)) - _eVirt.alpha;
      Eigen::MatrixXd amps = (_sss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += amps * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += amps.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd iajb = _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j)) - _eVirt.beta;
      Eigen::MatrixXd amps = (_sss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += amps * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += amps.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j)) - _eVirtab;
      Eigen::MatrixXd amps = (_oss * iajb).cwiseQuotient(denom);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += amps * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += amps.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  Eigen::MatrixXd rightHandSide = _Fia;
  Eigen::Map<Eigen::MatrixXd> rhs_a(rightHandSide.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> rhs_b(rightHandSide.data() + _alpha, _nv_b, _no_b);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
    rhs_a.noalias() -= zia_a * bij_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
    rhs_b.noalias() -= zia_b * bij_b.transpose();
  }
  Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nDim);
  rightHandSide.col(0).head(_alpha) += this->getJ2GContribution(_Zia.alpha.data(), nullptr, nullVector, true, 1);
  rightHandSide.col(0).tail(_beta) += this->getJ2GContribution(_Zia.beta.data(), nullptr, nullVector, true, -1);
  rightHandSide *= -1;

  return rightHandSide;
} /* this->calculateGroundStateLagrangeMultiplier() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity