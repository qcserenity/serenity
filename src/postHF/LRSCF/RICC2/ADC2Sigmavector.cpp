/**
 * @file ADC2Sigmavector.cpp
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
#include "postHF/LRSCF/RICC2/ADC2Controller.h"

/* Include Serenity Internal Headers */
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ADC2Controller<SCFMode>::ADC2Controller(std::shared_ptr<LRSCFController<SCFMode>> lrscf)
  : CC2Controller<SCFMode>(lrscf) {
}

template<>
void ADC2Controller<RESTRICTED>::calculateE() {
  _Eij = Eigen::MatrixXd::Zero(_no, _no);
  _Eab = Eigen::MatrixXd::Zero(_nv, _nv);

  _Eij.diagonal() = _e.segment(0, _no);
  _Eab.diagonal() = _e.segment(_no, _nv);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
    _Eij += jia.transpose() * yia;
    _Eab -= yia * jia.transpose();
  }

  // symmetrize
  _Eij = 0.5 * (_Eij + _Eij.transpose()).eval();
  _Eab = 0.5 * (_Eab + _Eab.transpose()).eval();
} /* this->calculateE() restricted */

template<>
void ADC2Controller<UNRESTRICTED>::calculateE() {
  auto& Jia = *_Jia;
  for_spin(_Eij, _Eab, _nv, _no, _Yia, Jia, _e) {
    _Eij_spin = Eigen::MatrixXd::Zero(_no_spin, _no_spin);
    _Eab_spin = Eigen::MatrixXd::Zero(_nv_spin, _nv_spin);

    _Eij_spin.diagonal() = _e_spin.segment(0, _no_spin);
    _Eab_spin.diagonal() = _e_spin.segment(_no_spin, _nv_spin);

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(Jia_spin.col(Q).data(), _nv_spin, _no_spin);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia_spin.col(Q).data(), _nv_spin, _no_spin);
      _Eij_spin += jia.transpose() * yia;
      _Eab_spin -= yia * jia.transpose();
    }

    // symmetrize
    _Eij_spin = 0.5 * (_Eij_spin + _Eij_spin.transpose()).eval();
    _Eab_spin = 0.5 * (_Eab_spin + _Eab_spin.transpose()).eval();
  };
} /* this->calculateE() unrestricted */

template<>
Eigen::VectorXd ADC2Controller<RESTRICTED>::getRightCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia, guessVector, false);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector;
  if (_settings.triplet) {
    sigmaVector = Eigen::VectorXd::Zero(_no * _nv);
  }
  else {
    sigmaVector = 2 * _Jia->operator*(_Jia->transpose() * guessVector);
  }
  Eigen::VectorXd Fia = sigmaVector;
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_no * _nv);

  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(Fia.data(), _nv, _no);

  sigma.noalias() += _Eab * guess - guess * _Eij;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    fia.noalias() -= jia * guess.transpose() * jia;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  Timings::takeTime("CC2 -         Amp-Contraction");
  _Zia.setZero();
  if (_settings.ltconv == 0) {
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = _settings.triplet ? 0 : i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
        Tia.segment(i * _nv, _nv).noalias() += tij * guessVector.segment(j * _nv, _nv);
        sigmaVector.segment(i * _nv, _nv).noalias() += 0.5 * tij * Fia.segment(j * _nv, _nv);
        _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
        if (i != j && !_settings.triplet) {
          Tia.segment(j * _nv, _nv).noalias() += tij.transpose() * guessVector.segment(i * _nv, _nv);
          sigmaVector.segment(j * _nv, _nv).noalias() += 0.5 * tij.transpose() * Fia.segment(i * _nv, _nv);
          _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
        }
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // MP2 amplitudes.
      Eigen::VectorXd param(_nv * _no);
      for (size_t i = 0, ia = 0; i < _no; ++i) {
        for (size_t a = _no; a < _no + _nv; ++a, ++ia) {
          param(ia) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
        }
      }

      Eigen::VectorXd X = (*_Jia).transpose() * (param.asDiagonal() * guessVector);
      Eigen::VectorXd Y = 0.5 * (*_Jia).transpose() * (param.asDiagonal() * Fia);

      Tia.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Jia) * X;
      sigmaVector.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Jia) * Y;

      // Right doubles amplitudes.
      param *= std::exp(0.5 * eigenvalue * _roots(iRoot));
      Eigen::MatrixXd V = ((*_Jia).transpose() * param.asDiagonal() * (*_Jia));
      Eigen::MatrixXd W = (_Xia.transpose() * param.asDiagonal() * (*_Jia));

      _Zia.noalias() -= param.asDiagonal() * _Xia * V;
      _Zia.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Jia) * W;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  if (!_settings.triplet) {
    sigmaVector.noalias() += 0.5 * 2 * _Jia->operator*(_Jia->transpose() * Tia);
  }
  Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij(_Jij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
    sigma.noalias() -= 0.5 * jia * tia.transpose() * jia;
    sigma.noalias() -= zia * jij;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  sigmaVector.noalias() += this->getJ2GContribution(_Zia.data(), _Jij->data(), guessVector, false);

  return sigmaVector;
} /* this->getRightCC2Sigma() restricted */

template<>
Eigen::VectorXd ADC2Controller<UNRESTRICTED>::getRightCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia.alpha, guessVector.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, guessVector.tail(_beta), false, -1);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * guessVector.head(_alpha) + _Jia->beta.transpose() * guessVector.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * XQ;
  Eigen::VectorXd Fia = sigmaVector;
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);

  Eigen::Map<Eigen::MatrixXd> guess_a(guessVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_a(Fia.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> guess_b(guessVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> fia_b(Fia.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

  sigma_a.noalias() += _Eab.alpha * guess_a - guess_a * _Eij.alpha;
  sigma_b.noalias() += _Eab.beta * guess_b - guess_b * _Eij.beta;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    fia_a.noalias() -= jia_a * guess_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    fia_b.noalias() -= jia_b * guess_b.transpose() * jia_b;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  Timings::takeTime("CC2 -         Amp-Contraction");
  _Zia.alpha.setZero();
  _Zia.beta.setZero();
  if (_settings.ltconv == 0) {
    // alpha
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 1);
        Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(j * _nv_a, _nv_a);
        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += 0.5 * tij * Fia.segment(j * _nv_a, _nv_a);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          Tia.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
          sigmaVector.segment(j * _nv_a, _nv_a).noalias() += 0.5 * tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
          _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }

    // beta
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
        Tia.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
        sigmaVector.segment(_alpha + i * _nv_b, _nv_b).noalias() += 0.5 * tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(_alpha + i * _nv_b, _nv_b);
          sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() +=
              0.5 * tij.transpose() * Fia.segment(_alpha + i * _nv_b, _nv_b);
          _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }

    // mixed
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);

        Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += 0.5 * tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);

        Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
        sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += 0.5 * tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // MP2 amplitudes.
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

      Eigen::VectorXd X_a = _Jia->alpha.transpose() * (param_a.asDiagonal() * guessVector.head(_alpha));
      Eigen::VectorXd X_b = _Jia->beta.transpose() * (param_b.asDiagonal() * guessVector.tail(_beta));
      Eigen::VectorXd Y_a = 0.5 * _Jia->alpha.transpose() * (param_a.asDiagonal() * Fia.head(_alpha));
      Eigen::VectorXd Y_b = 0.5 * _Jia->beta.transpose() * (param_b.asDiagonal() * Fia.tail(_beta));

      Tia.head(_alpha).noalias() -= param_a.asDiagonal() * _Jia->alpha * X_b;
      Tia.tail(_beta).noalias() -= param_b.asDiagonal() * _Jia->beta * X_a;
      sigmaVector.head(_alpha).noalias() -= param_a.asDiagonal() * _Jia->alpha * Y_b;
      sigmaVector.tail(_beta).noalias() -= param_b.asDiagonal() * _Jia->beta * Y_a;

      // Right doubles amplitudes.
      param_a *= std::exp(0.5 * eigenvalue * _roots(iRoot));
      param_b *= std::exp(0.5 * eigenvalue * _roots(iRoot));

      Eigen::MatrixXd V_a = (_Jia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha);
      Eigen::MatrixXd V_b = (_Jia->beta.transpose() * param_b.asDiagonal() * _Jia->beta);
      Eigen::MatrixXd W_a = (_Xia.alpha.transpose() * param_a.asDiagonal() * _Jia->alpha);
      Eigen::MatrixXd W_b = (_Xia.beta.transpose() * param_b.asDiagonal() * _Jia->beta);

      _Zia.alpha.noalias() -= param_a.asDiagonal() * _Xia.alpha * V_b;
      _Zia.beta.noalias() -= param_b.asDiagonal() * _Xia.beta * V_a;
      _Zia.alpha.noalias() -= param_a.asDiagonal() * _Jia->alpha * W_b;
      _Zia.beta.noalias() -= param_b.asDiagonal() * _Jia->beta * W_a;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
  sigmaVector.head(_alpha).noalias() += 0.5 * _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += 0.5 * _Jia->beta * XQ;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij_a(_Jij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
    sigma_a.noalias() -= 0.5 * jia_a * tia_a.transpose() * jia_a;
    sigma_a.noalias() -= zia_a * jij_a;

    Eigen::Map<Eigen::MatrixXd> jij_b(_Jij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
    sigma_b.noalias() -= 0.5 * jia_b * tia_b.transpose() * jia_b;
    sigma_b.noalias() -= zia_b * jij_b;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Zia.alpha.data(), _Jij->alpha.data(), guessVector.head(_alpha), false, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Zia.beta.data(), _Jij->beta.data(), guessVector.tail(_beta), false, -1);

  return sigmaVector;
} /* this->getRightCC2Sigma() unrestricted */

template class ADC2Controller<Options::SCF_MODES::RESTRICTED>;
template class ADC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity