/**
 * @file CC2Residual.cpp
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
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::calculateResidual() {
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);

  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> res(_residual.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);

  Eigen::VectorXd XQ = _Jia->transpose() * _singles;
  _Fia = 2 * _Jia->operator*(XQ);
  _residual = 2 * _Bia->operator*(XQ);

  res.noalias() += _e.segment(_no, _nv).asDiagonal() * sgl;
  res.noalias() -= sgl * _e.segment(0, _no).asDiagonal();

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    fia.noalias() -= jia * sgl.transpose() * jia;
  }

  _Yia.setZero();
  if (_settings.ltconv == 0) {
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j);
        _residual.segment(i * _nv, _nv).noalias() += tij * _Fia.segment(j * _nv, _nv);
        _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Jia->middleRows(j * _nv, _nv);
        if (i != j) {
          _residual.segment(j * _nv, _nv).noalias() += tij.transpose() * _Fia.segment(i * _nv, _nv);
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
          param(ia) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
        }
      }
      Eigen::VectorXd Y = (*_Bia).transpose() * (param.asDiagonal() * _Fia);
      _residual.noalias() -= param.asDiagonal() * (*_Bia) * Y;

      Eigen::MatrixXd V = ((*_Bia).transpose() * param.asDiagonal() * (*_Jia));
      _Yia.noalias() -= param.asDiagonal() * (*_Bia) * V;
    }
  }

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
    res.noalias() -= yia * bij;
  }

  _residual.noalias() += this->getJ2GContribution(_Yia.data(), _Bij->data(), _singles, false);
} /* this->calculateResidual() restricted */

template<>
void CC2Controller<UNRESTRICTED>::calculateResidual() {
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);

  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> res_a(_residual.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> res_b(_residual.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * _singles.head(_alpha) + _Jia->beta.transpose() * _singles.tail(_beta);
  _Fia.head(_alpha) += _Jia->alpha * XQ;
  _Fia.tail(_beta) += _Jia->beta * XQ;
  _residual.head(_alpha) += _Bia->alpha * XQ;
  _residual.tail(_beta) += _Bia->beta * XQ;

  res_a.noalias() += _e.alpha.segment(_no_a, _nv_a).asDiagonal() * sgl_a;
  res_a.noalias() -= sgl_a * _e.alpha.segment(0, _no_a).asDiagonal();

  res_b.noalias() += _e.beta.segment(_no_b, _nv_b).asDiagonal() * sgl_b;
  res_b.noalias() -= sgl_b * _e.beta.segment(0, _no_b).asDiagonal();

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    fia_a.noalias() -= jia_a * sgl_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    fia_b.noalias() -= jia_b * sgl_b.transpose() * jia_b;
  }
  _Yia.alpha.setZero();
  _Yia.beta.setZero();
  if (_settings.ltconv == 0) {
    // alpha
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
        _residual.segment(i * _nv_a, _nv_a).noalias() += tij * _Fia.segment(j * _nv_a, _nv_a);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _residual.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }

    // beta
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
        _residual.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _residual.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * _Fia.segment(_alpha + i * _nv_b, _nv_b);
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }

    // mixed
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
        _residual.segment(i * _nv_a, _nv_a).noalias() += tij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);

        _residual.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      Eigen::VectorXd param_a(_alpha);
      for (size_t i = 0, ia = 0; i < _no_a; ++i) {
        for (size_t a = _no_a; a < _no_a + _nv_a; ++a, ++ia) {
          param_a(ia) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e.alpha(a) - _e.alpha(i)) * _roots(iRoot));
        }
      }
      Eigen::VectorXd param_b(_beta);
      for (size_t j = 0, jb = 0; j < _no_b; ++j) {
        for (size_t b = _no_b; b < _no_b + _nv_b; ++b, ++jb) {
          param_b(jb) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e.beta(b) - _e.beta(j)) * _roots(iRoot));
        }
      }
      Eigen::VectorXd Y_a = _Bia->alpha.transpose() * (param_a.asDiagonal() * _Fia.head(_alpha));
      Eigen::VectorXd Y_b = _Bia->beta.transpose() * (param_b.asDiagonal() * _Fia.tail(_beta));

      _residual.head(_alpha).noalias() -= param_a.asDiagonal() * _Bia->alpha * Y_b;
      _residual.tail(_beta).noalias() -= param_b.asDiagonal() * _Bia->beta * Y_a;

      Eigen::MatrixXd V_a = (_Bia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha);
      Eigen::MatrixXd V_b = (_Bia->beta.transpose() * param_b.asDiagonal() * _Jia->beta);
      _Yia.alpha.noalias() -= param_a.asDiagonal() * _Bia->alpha * V_b;
      _Yia.beta.noalias() -= param_b.asDiagonal() * _Bia->beta * V_a;
    }
  }

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);
    res_a.noalias() -= yia_a * bij_a;

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);
    res_b.noalias() -= yia_b * bij_b;
  }

  _residual.head(_alpha).noalias() +=
      this->getJ2GContribution(_Yia.alpha.data(), _Bij->alpha.data(), _singles.head(_alpha), false, 1);
  _residual.tail(_beta).noalias() +=
      this->getJ2GContribution(_Yia.beta.data(), _Bij->beta.data(), _singles.tail(_beta), false, -1);
} /* this->calculateResidual() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity