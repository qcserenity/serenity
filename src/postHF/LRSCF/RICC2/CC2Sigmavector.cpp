/**
 * @file CC2Sigmavector.cpp
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
void CC2Controller<RESTRICTED>::calculateE() {
  _Eij = Eigen::MatrixXd::Zero(_no, _no);
  _Eab = Eigen::MatrixXd::Zero(_nv, _nv);

  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the orbital energy differences: (e_a - e_i) b_{ai}.
  if (_xwfModel != Options::LR_METHOD::CISD) {
    _Eij.diagonal() = _e.segment(0, _no);
    _Eab.diagonal() = _e.segment(_no, _nv);
  }

  Eigen::VectorXd XQ = _Jia->transpose() * _singles;
  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
    // ^Fij and ^Fab
    _Eij += 2 * bij * XQ(Q);

    _Eij -= jia.transpose() * sgl * bij;

    _Eij += jia.transpose() * yia;
    _Eab -= yia * jia.transpose();
  }

  _Eab += this->transformedABFock(_singles, XQ, _Yia);
} /* this->calculateE() restricted */

template<>
void CC2Controller<UNRESTRICTED>::calculateE() {
  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);

  _Eij.alpha = Eigen::MatrixXd::Zero(_no_a, _no_a);
  _Eij.beta = Eigen::MatrixXd::Zero(_no_b, _no_b);
  _Eab.alpha = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
  _Eab.beta = Eigen::MatrixXd::Zero(_nv_b, _nv_b);

  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the orbital energy differences: (e_a - e_i) b_{ai}.
  if (_xwfModel != Options::LR_METHOD::CISD) {
    _Eij.alpha.diagonal() = _e.alpha.segment(0, _no_a);
    _Eij.beta.diagonal() = _e.beta.segment(0, _no_b);
    _Eab.alpha.diagonal() = _e.alpha.segment(_no_a, _nv_a);
    _Eab.beta.diagonal() = _e.beta.segment(_no_b, _nv_b);
  }

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * _singles.head(_alpha) + _Jia->beta.transpose() * _singles.tail(_beta);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);
    // ^Fij and ^Fab
    _Eij.alpha += bij_a * XQ(Q);

    _Eij.alpha -= jia_a.transpose() * sgl_a * bij_a;

    _Eij.alpha += jia_a.transpose() * yia_a;
    _Eab.alpha -= yia_a * jia_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);
    // ^Fij and ^Fab
    _Eij.beta += bij_b * XQ(Q);

    _Eij.beta -= jia_b.transpose() * sgl_b * bij_b;

    _Eij.beta += jia_b.transpose() * yia_b;
    _Eab.beta -= yia_b * jia_b.transpose();
  }

  _Eab.alpha += this->transformedABFock(_singles.head(_alpha), XQ, _Yia.alpha, 1);
  _Eab.beta += this->transformedABFock(_singles.tail(_beta), XQ, _Yia.beta, -1);
} /* this->calculateE() unrestricted */

template<>
Eigen::VectorXd CC2Controller<RESTRICTED>::getRightCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia, guessVector, false);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Fia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Jia->transpose() * guessVector;
  if (!_settings.triplet) {
    sigmaVector = 2 * _Bia->operator*(XQ);
    Fia = 2 * _Jia->operator*(XQ);
  }
  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the Coulomb term: (ai|jb) b_{bj}.
  if (_xwfModel == Options::LR_METHOD::CISD) {
    sigmaVector.setZero();
  }

  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
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
        sigmaVector.segment(i * _nv, _nv).noalias() += tij * Fia.segment(j * _nv, _nv);
        sigmaVector.segment(i * _nv, _nv).noalias() += rij * _Fia.segment(j * _nv, _nv);
        _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
        if (i != j && !_settings.triplet) {
          sigmaVector.segment(j * _nv, _nv).noalias() += tij.transpose() * Fia.segment(i * _nv, _nv);
          sigmaVector.segment(j * _nv, _nv).noalias() += rij.transpose() * _Fia.segment(i * _nv, _nv);
          _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
        }
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // CC2 amplitudes.
      Eigen::VectorXd param(_nv * _no);
      for (size_t i = 0, ia = 0; i < _no; ++i) {
        for (size_t a = _no; a < _no + _nv; ++a, ++ia) {
          param(ia) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
        }
      }
      Eigen::VectorXd X = (*_Bia).transpose() * (param.asDiagonal() * Fia);
      sigmaVector.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Bia) * X;

      // Right doubles amplitudes.
      param *= std::exp(0.5 * eigenvalue * _roots(iRoot));

      Eigen::VectorXd Y = (*_Bia).transpose() * (param.asDiagonal() * _Fia);
      Eigen::VectorXd Z = _Xia.transpose() * (param.asDiagonal() * _Fia);

      sigmaVector.noalias() -= param.asDiagonal() * _Xia * Y;
      sigmaVector.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Bia) * Z;

      Eigen::MatrixXd V = ((*_Bia).transpose() * param.asDiagonal() * (*_Jia));
      Eigen::MatrixXd W = (_Xia.transpose() * param.asDiagonal() * (*_Jia));

      _Zia.noalias() -= param.asDiagonal() * _Xia * V;
      _Zia.noalias() -= (_settings.triplet ? -1.0 : 1.0) * param.asDiagonal() * (*_Bia) * W;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
    sigma.noalias() -= zia * bij;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the exchange term: (ab|ji) b_{bj}.
  // To remove it, just pass a nullptr for the (Q|ji) integrals.
  double* ijpointer = (_xwfModel == Options::LR_METHOD::CISD) ? nullptr : _Bij->data();
  sigmaVector += this->getJ2GContribution(_Zia.data(), ijpointer, guessVector, false);

  return sigmaVector;
} /* this->getRightCC2Sigma() restricted */

template<>
Eigen::VectorXd CC2Controller<UNRESTRICTED>::getRightCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia.alpha, guessVector.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, guessVector.tail(_beta), false, -1);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Fia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * guessVector.head(_alpha) + _Jia->beta.transpose() * guessVector.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Bia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Bia->beta * XQ;
  Fia.head(_alpha).noalias() += _Jia->alpha * XQ;
  Fia.tail(_beta).noalias() += _Jia->beta * XQ;

  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the Coulomb term: (ai|jb) b_{bj}.
  if (_xwfModel == Options::LR_METHOD::CISD) {
    sigmaVector.setZero();
  }

  Eigen::Map<Eigen::MatrixXd> guess_a(guessVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_a(Fia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> guess_b(guessVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> fia_b(Fia.data() + _alpha, _nv_b, _no_b);

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
        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += tij * Fia.segment(j * _nv_a, _nv_a);
        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += rij * _Fia.segment(j * _nv_a, _nv_a);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          sigmaVector.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
          sigmaVector.segment(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
          _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }

    // beta
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
        sigmaVector.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
        sigmaVector.segment(_alpha + i * _nv_b, _nv_b).noalias() += rij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * Fia.segment(_alpha + i * _nv_b, _nv_b);
          sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() +=
              rij.transpose() * _Fia.segment(_alpha + i * _nv_b, _nv_b);
          _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }

    // mixed
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);

        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
        sigmaVector.segment(i * _nv_a, _nv_a).noalias() += rij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);

        sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
        sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += rij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // CC2 amplitudes.
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
      Eigen::VectorXd X_a = _Bia->alpha.transpose() * (param_a.asDiagonal() * Fia.head(_alpha));
      Eigen::VectorXd X_b = _Bia->beta.transpose() * (param_b.asDiagonal() * Fia.tail(_beta));
      sigmaVector.head(_alpha).noalias() -= param_a.asDiagonal() * _Bia->alpha * X_b;
      sigmaVector.tail(_beta).noalias() -= param_b.asDiagonal() * _Bia->beta * X_a;

      // Right doubles amplitudes.
      param_a *= std::exp(0.5 * eigenvalue * _roots(iRoot));
      param_b *= std::exp(0.5 * eigenvalue * _roots(iRoot));

      Eigen::VectorXd Y_a = _Bia->alpha.transpose() * (param_a.asDiagonal() * _Fia.head(_alpha));
      Eigen::VectorXd Y_b = _Bia->beta.transpose() * (param_b.asDiagonal() * _Fia.tail(_beta));
      Eigen::VectorXd Z_a = _Xia.alpha.transpose() * (param_a.asDiagonal() * _Fia.head(_alpha));
      Eigen::VectorXd Z_b = _Xia.beta.transpose() * (param_b.asDiagonal() * _Fia.tail(_beta));

      sigmaVector.head(_alpha).noalias() -= param_a.asDiagonal() * _Xia.alpha * Y_b;
      sigmaVector.tail(_beta).noalias() -= param_b.asDiagonal() * _Xia.beta * Y_a;
      sigmaVector.head(_alpha).noalias() -= param_a.asDiagonal() * _Bia->alpha * Z_b;
      sigmaVector.tail(_beta).noalias() -= param_b.asDiagonal() * _Bia->beta * Z_a;

      Eigen::MatrixXd V_a = (_Bia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha);
      Eigen::MatrixXd V_b = (_Bia->beta.transpose() * param_b.asDiagonal() * _Jia->beta);
      Eigen::MatrixXd W_a = (_Xia.alpha.transpose() * param_a.asDiagonal() * _Jia->alpha);
      Eigen::MatrixXd W_b = (_Xia.beta.transpose() * param_b.asDiagonal() * _Jia->beta);

      _Zia.alpha.noalias() -= param_a.asDiagonal() * _Xia.alpha * V_b;
      _Zia.beta.noalias() -= param_b.asDiagonal() * _Xia.beta * V_a;
      _Zia.alpha.noalias() -= param_a.asDiagonal() * _Bia->alpha * W_b;
      _Zia.beta.noalias() -= param_b.asDiagonal() * _Bia->beta * W_a;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
    sigma_a.noalias() -= zia_a * bij_a;

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
    sigma_b.noalias() -= zia_b * bij_b;
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  // If only the CIS(D) correction is requested, all sigma vector terms
  // corresponding to the CIS matrix need to be removed from this sigma build.
  // This one corresponds to the exchange term: (ab|ji) b_{bj}.
  // To remove it, just pass a nullptr for the (Q|ji) integrals.
  double* ijpointer_a = (_xwfModel == Options::LR_METHOD::CISD) ? nullptr : _Bij->alpha.data();
  double* ijpointer_b = (_xwfModel == Options::LR_METHOD::CISD) ? nullptr : _Bij->beta.data();
  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Zia.alpha.data(), ijpointer_a, guessVector.head(_alpha), false, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Zia.beta.data(), ijpointer_b, guessVector.tail(_beta), false, -1);

  return sigmaVector;
} /* this->getRightCC2Sigma() unrestricted */

template<>
Eigen::VectorXd CC2Controller<RESTRICTED>::getLeftCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  _leftSingles = guessVector;
  this->performTransformation(_Zia, guessVector, true);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector = 2 * _Jia->operator*(_Bia->transpose() * guessVector);

  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);

  sigma.noalias() += _Eab.transpose() * guess - guess * _Eij.transpose();
  Timings::timeTaken("CC2 -           Q-Contraction");

  Timings::takeTime("CC2 -         Amp-Contraction");
  _Xia.setZero();
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_no * _nv);
  if (_settings.ltconv == 0) {
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
        Tia.segment(i * _nv, _nv).noalias() += tij * guessVector.segment(j * _nv, _nv);
        _Xia.middleRows(i * _nv, _nv).noalias() += lij * _Bia->middleRows(j * _nv, _nv);
        if (i != j) {
          Tia.segment(j * _nv, _nv).noalias() += tij.transpose() * guessVector.segment(i * _nv, _nv);
          _Xia.middleRows(j * _nv, _nv).noalias() += lij.transpose() * _Bia->middleRows(i * _nv, _nv);
        }
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // CC2 amplitudes.
      Eigen::VectorXd param(_nv * _no);
      for (unsigned i = 0, ia = 0; i < _no; ++i) {
        for (unsigned a = _no; a < _no + _nv; ++a, ++ia) {
          param(ia) = std::sqrt(_oss * _weights(iRoot)) * std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
        }
      }

      Eigen::VectorXd X = (*_Bia).transpose() * (param.asDiagonal() * guessVector);
      Tia.noalias() -= param.asDiagonal() * (*_Bia) * X;

      // Left doubles amplitudes.
      param *= std::exp(0.5 * eigenvalue * _roots(iRoot));
      Eigen::MatrixXd V = ((*_Jia).transpose() * param.asDiagonal() * (*_Bia));
      Eigen::MatrixXd W = (_Zia.transpose() * param.asDiagonal() * (*_Bia));
      Eigen::MatrixXd Y = _Fia.transpose() * param.asDiagonal() * (*_Bia);
      Eigen::MatrixXd Z = _leftSingles.transpose() * param.asDiagonal() * (*_Bia);

      _Xia.noalias() -= param.asDiagonal() * _Zia * V;
      _Xia.noalias() -= param.asDiagonal() * (*_Jia) * W;
      _Xia.noalias() -= param.asDiagonal() * _leftSingles * Y;
      _Xia.noalias() -= param.asDiagonal() * _Fia * Z;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  sigmaVector.noalias() += 2 * _Jia->operator*(_Jia->transpose() * Tia);
  Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> xia(_Xia.col(Q).data(), _nv, _no);
    sigma.noalias() -= jia * tia.transpose() * jia;
    sigma.noalias() -= xia * bij.transpose();
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  sigmaVector += this->getJ2GContribution(_Xia.data(), _Bij->data(), guessVector, true);

  return sigmaVector;
} /* this->getLeftCC2Sigma() restricted */

template<>
Eigen::VectorXd CC2Controller<UNRESTRICTED>::getLeftCC2Sigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  _leftSingles = guessVector;
  this->performTransformation(_Zia.alpha, guessVector.head(_alpha), true, 1);
  this->performTransformation(_Zia.beta, guessVector.tail(_beta), true, -1);

  Timings::takeTime("CC2 -           Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Bia->alpha.transpose() * guessVector.head(_alpha) + _Bia->beta.transpose() * guessVector.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * XQ;

  Eigen::Map<Eigen::MatrixXd> guess_a(guessVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> guess_b(guessVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

  sigma_a.noalias() += _Eab.alpha.transpose() * guess_a - guess_a * _Eij.alpha.transpose();
  sigma_b.noalias() += _Eab.beta.transpose() * guess_b - guess_b * _Eij.beta.transpose();
  Timings::timeTaken("CC2 -           Q-Contraction");

  Timings::takeTime("CC2 -         Amp-Contraction");
  _Xia.alpha.setZero();
  _Xia.beta.setZero();
  if (_settings.ltconv == 0) {
    // alpha
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
        Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(j * _nv_a, _nv_a);
        _Xia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += lij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          Tia.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
          _Xia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += lij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }

    // beta
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
        Tia.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
        _Xia.beta.middleRows(i * _nv_b, _nv_b).noalias() += lij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(_alpha + i * _nv_b, _nv_b);
          _Xia.beta.middleRows(j * _nv_b, _nv_b).noalias() += lij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }

    // mixed
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);

        Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
        _Xia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += lij * _Bia->beta.middleRows(j * _nv_b, _nv_b);

        Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
        _Xia.beta.middleRows(j * _nv_b, _nv_b).noalias() += lij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  else {
    for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
      // CC2 amplitudes.
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

      Eigen::VectorXd X_a = _Bia->alpha.transpose() * (param_a.asDiagonal() * guessVector.head(_alpha));
      Eigen::VectorXd X_b = _Bia->beta.transpose() * (param_b.asDiagonal() * guessVector.tail(_beta));
      Tia.head(_alpha).noalias() -= param_a.asDiagonal() * _Bia->alpha * X_b;
      Tia.tail(_beta).noalias() -= param_b.asDiagonal() * _Bia->beta * X_a;

      // Left doubles amplitudes.
      param_a *= std::exp(0.5 * eigenvalue * _roots(iRoot));
      param_b *= std::exp(0.5 * eigenvalue * _roots(iRoot));

      Eigen::MatrixXd V_a = (_Jia->alpha.transpose() * param_a.asDiagonal() * _Bia->alpha);
      Eigen::MatrixXd V_b = (_Jia->beta.transpose() * param_b.asDiagonal() * _Bia->beta);
      Eigen::MatrixXd W_a = (_Zia.alpha.transpose() * param_a.asDiagonal() * _Bia->alpha);
      Eigen::MatrixXd W_b = (_Zia.beta.transpose() * param_b.asDiagonal() * _Bia->beta);
      Eigen::MatrixXd Y_a = _Fia.head(_alpha).transpose() * param_a.asDiagonal() * _Bia->alpha;
      Eigen::MatrixXd Y_b = _Fia.tail(_beta).transpose() * param_b.asDiagonal() * _Bia->beta;
      Eigen::MatrixXd Z_a = _leftSingles.head(_alpha).transpose() * param_a.asDiagonal() * _Bia->alpha;
      Eigen::MatrixXd Z_b = _leftSingles.tail(_beta).transpose() * param_b.asDiagonal() * _Bia->beta;

      _Xia.alpha.noalias() -= param_a.asDiagonal() * _Zia.alpha * V_b;
      _Xia.beta.noalias() -= param_b.asDiagonal() * _Zia.beta * V_a;
      _Xia.alpha.noalias() -= param_a.asDiagonal() * _Jia->alpha * W_b;
      _Xia.beta.noalias() -= param_b.asDiagonal() * _Jia->beta * W_a;
      _Xia.alpha.noalias() -= param_a.asDiagonal() * _leftSingles.head(_alpha) * Y_b;
      _Xia.beta.noalias() -= param_b.asDiagonal() * _leftSingles.tail(_beta) * Y_a;
      _Xia.alpha.noalias() -= param_a.asDiagonal() * _Fia.head(_alpha) * Z_b;
      _Xia.beta.noalias() -= param_b.asDiagonal() * _Fia.tail(_beta) * Z_a;
    }
  }
  Timings::timeTaken("CC2 -         Amp-Contraction");

  Timings::takeTime("CC2 -           Q-Contraction");
  XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * XQ;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.col(Q).data(), _nv_a, _no_a);
    sigma_a.noalias() -= jia_a * tia_a.transpose() * jia_a;
    sigma_a.noalias() -= xia_a * bij_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.col(Q).data(), _nv_b, _no_b);
    sigma_b.noalias() -= jia_b * tia_b.transpose() * jia_b;
    sigma_b.noalias() -= xia_b * bij_b.transpose();
  }
  Timings::timeTaken("CC2 -           Q-Contraction");

  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Xia.alpha.data(), _Bij->alpha.data(), guessVector.head(_alpha), true, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Xia.beta.data(), _Bij->beta.data(), guessVector.tail(_beta), true, -1);

  return sigmaVector;
} /* this->getLeftCC2Sigma() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity