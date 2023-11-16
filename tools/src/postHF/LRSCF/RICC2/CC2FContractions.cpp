/**
 * @file CC2FContractions.cpp
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

namespace Serenity {

template<>
Eigen::VectorXd CC2Controller<RESTRICTED>::getFSingles(Eigen::Ref<Eigen::VectorXd> amp) {
  // Maps.
  Eigen::Map<Eigen::MatrixXd> t(amp.data(), _nv, _no);

  Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);

  Eigen::MatrixXd Eij = +1.0 * fia.transpose() * t;
  Eigen::MatrixXd Eab = -1.0 * t * fia.transpose();

  Eigen::VectorXd XQ = _Jia->transpose() * amp;

  // 01 contribution.
  Eigen::VectorXd sigmaVector = 2 * (*_Jia) * XQ;
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);

    Eij.noalias() += 2 * bij * XQ(Q);
    Eij.noalias() -= jia.transpose() * t * bij;

    // 02 contribution.
    sigma.noalias() -= jia * t.transpose() * jia;
  }

  Eab += this->transformedABFock(amp, XQ, _Zia);

  // 03 and 04 contribution.
  sigma += Eab.transpose() * gsl - gsl * Eij.transpose();

  this->performTransformation(_Xia, amp, false);
  this->performTransformation(_Zia, _gsLagrange, true);

  // Amplitude contractions for F and J contributions.
  _Yia = Eigen::MatrixXd::Zero(_nv * _no, _nx);                // hat Y
  Eigen::MatrixXd Zia = Eigen::MatrixXd::Zero(_nv * _no, _nx); // bar Y
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
      Zia.middleRows(i * _nv, _nv).noalias() += tij * _Xia.middleRows(j * _nv, _nv);
      _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Bia->middleRows(j * _nv, _nv);
      if (i != j) {
        Zia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Xia.middleRows(i * _nv, _nv);
        _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Bia->middleRows(i * _nv, _nv);
      }
    }
  }

  // J1 contribution.
  sigmaVector += 2 * (*_Jia) * (_Xia.transpose() * _gsLagrange);

  // F and J contrations.
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> zia(Zia.col(Q).data(), _nv, _no);

    Eigen::MatrixXd _bab = jia * t.transpose();
    Eigen::MatrixXd _bij = t.transpose() * jia;

    // F2 contribution.
    sigma.noalias() -= zia * bij.transpose();

    // F3 contribution.
    sigma.noalias() -= _bab * yia;

    // F4 contribution.
    sigma.noalias() -= yia * _bij;

    // J3 contribution.
    sigma.noalias() += _bab * gsl * bij.transpose();

    // Prepare J2 contribution.
    zia.noalias() -= gsl * _bij;
  }

  // F1 and J2 contribution.
  Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nv * _no);
  sigmaVector += this->getJ2GContribution(Zia.data(), nullptr, nullVector, true);

  return sigmaVector;
}

template<>
Eigen::VectorXd CC2Controller<RESTRICTED>::getFDoubles(Eigen::Ref<Eigen::VectorXd> amp, Eigen::Ref<Eigen::VectorXd> dip,
                                                       double frequency) {
  // Dipole Integrals
  unsigned mo = _no + _nv;
  Eigen::Map<Eigen::MatrixXd> D(dip.data(), mo, mo);
  SPMatrix<RESTRICTED> dip_oo = D.topLeftCorner(_no, _no);
  SPMatrix<RESTRICTED> dip_vv = D.bottomRightCorner(_nv, _nv);

  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nv * _no);
  Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

  this->performTransformation(_Xia, amp, false);
  _Zia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, frequency) + this->getXiDoubles(i, j, frequency, dip_vv);
      _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
      Tia.segment(i * _nv, _nv).noalias() += rij * _gsLagrange.segment(j * _nv, _nv);
      if (i != j) {
        _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
        Tia.segment(j * _nv, _nv).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv, _nv);
      }
    }
  }

  this->reorderTensor(*_Jia, _nv, _no);
  if (_Bia != _Jia) {
    this->reorderTensor(*_Bia, _nv, _no);
  }
  this->reorderTensor(_Zia, _nv, _no);
  this->reorderTensor(Tia, _nv, _no);
  this->reorderTensor(_gsLagrange, _nv, _no);

  for (size_t a = 0; a < _nv; ++a) {
    for (size_t b = a; b < _nv; ++b) {
      Eigen::MatrixXd rab = this->getXiDoublesV(a, b, frequency, dip_oo);
      _Zia.middleRows(a * _no, _no).noalias() += rab * _Jia->middleRows(b * _no, _no);
      Tia.segment(a * _no, _no).noalias() += rab * _gsLagrange.segment(b * _no, _no);
      if (a != b) {
        _Zia.middleRows(b * _no, _no).noalias() += rab.transpose() * _Jia->middleRows(a * _no, _no);
        Tia.segment(b * _no, _no).noalias() += rab.transpose() * _gsLagrange.segment(a * _no, _no);
      }
    }
  }

  this->reorderTensor(*_Jia, _no, _nv);
  if (_Bia != _Jia) {
    this->reorderTensor(*_Bia, _no, _nv);
  }
  this->reorderTensor(_Zia, _no, _nv);
  this->reorderTensor(Tia, _no, _nv);
  this->reorderTensor(_gsLagrange, _no, _nv);

  Eigen::VectorXd sigmaVector = 2 * _Jia->operator*(_Jia->transpose() * Tia);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);

  // GH contribution
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    sigma.noalias() -= jia * tia.transpose() * jia;
  }

  return sigmaVector;
}

template<>
std::pair<std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>>
CC2Controller<RESTRICTED>::getFContractions(std::vector<Eigen::MatrixXd>& electricAmps,
                                            std::vector<Eigen::MatrixXd>& magneticAmps, std::vector<double> frequencies,
                                            Eigen::Ref<Eigen::MatrixXd> electricDips,
                                            Eigen::Ref<Eigen::MatrixXd> magneticDips, Options::GAUGE gauge) {
  unsigned nFreqs = frequencies.size();
  std::pair<std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>> dotProducts;

  dotProducts.first = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());
  dotProducts.second = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());

  double gaugeFactorForORD = (gauge == Options::GAUGE::VELOCITY) ? 1 : -1;

  auto dressIntegrals = [&](Eigen::Ref<Eigen::VectorXd> amp) {
    Timings::takeTime("CC2 -              Dress Ints");
    double dot = 0;
    Eigen::Map<Eigen::MatrixXd> t(amp.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
    Eigen::MatrixXd tg = t.transpose() * gsl;
    Eigen::MatrixXd gt = gsl * t.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
      Eigen::MatrixXd tmp = jia * tg + gt * jia;
      dot += zia.cwiseProduct(tmp).sum();
    }
    Timings::timeTaken("CC2 -              Dress Ints");
    return dot;
  };

  Timings::takeTime("CC2 -         F-Matrix Contr.");
  for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
    double frequency = frequencies[iFreq];
    unsigned p = iFreq;
    unsigned m = (frequencies[iFreq] != 0) ? iFreq + nFreqs : iFreq;
    double F = (frequencies[iFreq] == 0) ? 2.0 : 1.0;
    for (unsigned I = 0; I < 3; ++I) {
      Eigen::Ref<Eigen::VectorXd> epI = electricAmps[p].col(I);
      Eigen::Ref<Eigen::VectorXd> emI = electricAmps[m].col(I);

      Eigen::Ref<Eigen::VectorXd> mpI = magneticAmps[p].col(I);
      Eigen::Ref<Eigen::VectorXd> mmI = magneticAmps[m].col(I);

      /***********/
      /* Singles */
      /***********/
      Eigen::VectorXd eps = this->getFSingles(epI);
      // <r+|r->
      dotProducts.first[iFreq].row(I) += F * eps.transpose() * electricAmps[m];
      // <r+|m->
      dotProducts.second[iFreq].row(I) += gaugeFactorForORD * F * eps.transpose() * magneticAmps[m];

      if (frequency != 0) {
        Eigen::VectorXd ems = this->getFSingles(emI);
        // <r-|r+>
        dotProducts.first[iFreq].row(I) += ems.transpose() * electricAmps[p];
        // <r-|m+>
        dotProducts.second[iFreq].row(I) += ems.transpose() * magneticAmps[p];
      }

      /***********/
      /* Doubles */
      /***********/
      // Electric amplitudes plus frequency.
      auto epd = this->getFDoubles(epI, electricDips.col(I), frequency);
      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> emJ = electricAmps[m].col(J);
        double dot = epd.dot(emJ) - dressIntegrals(emJ);
        // <r+|r->
        dotProducts.first[iFreq](I, J) += F * dot;
        dotProducts.first[iFreq](J, I) += F * dot;
      }

      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> mmJ = magneticAmps[m].col(J);
        double dot = epd.dot(mmJ) - dressIntegrals(mmJ);
        // <r+|m->
        dotProducts.second[iFreq](I, J) += gaugeFactorForORD * F * dot;
      }

      // Magnetic amplitudes plus frequency.
      auto mpd = this->getFDoubles(mpI, magneticDips.col(I), frequency);
      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> emJ = electricAmps[m].col(J);
        double dot = mpd.dot(emJ) - dressIntegrals(emJ);
        // <r-|m+>
        dotProducts.second[iFreq](J, I) += F * dot;
      }

      if (frequency != 0) {
        // Electric amplitudes minus frequency.
        auto emd = this->getFDoubles(emI, electricDips.col(I), -frequency);
        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> epJ = electricAmps[p].col(J);
          double dot = emd.dot(epJ) - dressIntegrals(epJ);
          // <r-|r+>
          dotProducts.first[iFreq](I, J) += dot;
          dotProducts.first[iFreq](J, I) += dot;
        }

        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> mpJ = magneticAmps[p].col(J);
          double dot = emd.dot(mpJ) - dressIntegrals(mpJ);
          // <r-|m+>
          dotProducts.second[iFreq](I, J) += dot;
        }

        // Magnetic amplitudes minus frequency.
        auto mmd = this->getFDoubles(mmI, magneticDips.col(I), -frequency);
        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> epJ = electricAmps[p].col(J);
          double dot = mmd.dot(epJ) - dressIntegrals(epJ);
          // <r+|m->
          dotProducts.second[iFreq](J, I) += gaugeFactorForORD * dot;
        }
      }
    }
  } /* for iFreq */
  Timings::timeTaken("CC2 -         F-Matrix Contr.");

  return dotProducts;
} /* this->calcFContractions() restricted */

template<>
Eigen::VectorXd CC2Controller<UNRESTRICTED>::getFSingles(Eigen::Ref<Eigen::VectorXd> amp) {
  // Maps.
  Eigen::Map<Eigen::MatrixXd> t_a(amp.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> t_b(amp.data() + _alpha, _nv_b, _no_b);

  Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);

  Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);

  Eigen::MatrixXd Eij_a = +1.0 * fia_a.transpose() * t_a;
  Eigen::MatrixXd Eij_b = +1.0 * fia_b.transpose() * t_b;
  Eigen::MatrixXd Eab_a = -1.0 * t_a * fia_a.transpose();
  Eigen::MatrixXd Eab_b = -1.0 * t_b * fia_b.transpose();

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * amp.head(_alpha) + _Jia->beta.transpose() * amp.tail(_beta);

  // 01 contribution.
  Eigen::VectorXd sigmaVector(_nDim);
  sigmaVector.head(_alpha) = _Jia->alpha * XQ;
  sigmaVector.tail(_beta) = _Jia->beta * XQ;

  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);

    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);

    Eij_a.noalias() += bij_a * XQ(Q);
    Eij_b.noalias() += bij_b * XQ(Q);
    Eij_a.noalias() -= jia_a.transpose() * t_a * bij_a;
    Eij_b.noalias() -= jia_b.transpose() * t_b * bij_b;

    // 02 contribution.
    sigma_a.noalias() -= jia_a * t_a.transpose() * jia_a;
    sigma_b.noalias() -= jia_b * t_b.transpose() * jia_b;
  }

  Eab_a += this->transformedABFock(amp.head(_alpha), XQ, _Zia.alpha, 1);
  Eab_b += this->transformedABFock(amp.tail(_beta), XQ, _Zia.beta, -1);

  // 03 and 04 contribution.
  sigma_a += Eab_a.transpose() * gsl_a - gsl_a * Eij_a.transpose();
  sigma_b += Eab_b.transpose() * gsl_b - gsl_b * Eij_b.transpose();

  this->performTransformation(_Xia.alpha, amp.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, amp.tail(_beta), false, -1);
  this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
  this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);

  // Amplitude contractions for F and J contributions.
  _Yia.alpha = Eigen::MatrixXd::Zero(_nv_a * _no_a, _nx);            // hat Y
  _Yia.beta = Eigen::MatrixXd::Zero(_nv_b * _no_b, _nx);             // hat Y
  Eigen::MatrixXd Zia_a = Eigen::MatrixXd::Zero(_nv_a * _no_a, _nx); // bar Y
  Eigen::MatrixXd Zia_b = Eigen::MatrixXd::Zero(_nv_b * _no_b, _nx); // bar Y

  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
      Zia_a.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.alpha.middleRows(j * _nv_a, _nv_a);
      _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        Zia_a.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
        _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
      Zia_b.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
      _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        Zia_b.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.beta.middleRows(i * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
      Zia_a.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
      _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      Zia_b.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
      _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  // J1 contribution.
  Eigen::VectorXd YQ = _Xia.alpha.transpose() * _gsLagrange.head(_alpha) + _Xia.beta.transpose() * _gsLagrange.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * YQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * YQ;

  // F and J contrations.
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(Zia_a.col(Q).data(), _nv_a, _no_a);

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(Zia_b.col(Q).data(), _nv_b, _no_b);

    Eigen::MatrixXd _bab_a = jia_a * t_a.transpose();
    Eigen::MatrixXd _bab_b = jia_b * t_b.transpose();
    Eigen::MatrixXd _bij_a = t_a.transpose() * jia_a;
    Eigen::MatrixXd _bij_b = t_b.transpose() * jia_b;

    // F2 contribution.
    sigma_a.noalias() -= zia_a * bij_a.transpose();
    sigma_b.noalias() -= zia_b * bij_b.transpose();

    // F3 contribution.
    sigma_a.noalias() -= _bab_a * yia_a;
    sigma_b.noalias() -= _bab_b * yia_b;

    // F4 contribution.
    sigma_a.noalias() -= yia_a * _bij_a;
    sigma_b.noalias() -= yia_b * _bij_b;

    // J3 contribution.
    sigma_a.noalias() += _bab_a * gsl_a * bij_a.transpose();
    sigma_b.noalias() += _bab_b * gsl_b * bij_b.transpose();

    // Prepare J2 contribution.
    zia_a.noalias() -= gsl_a * _bij_a;
    zia_b.noalias() -= gsl_b * _bij_b;
  }

  // F1 and J2 contribution.
  Eigen::VectorXd nullVector_a = Eigen::VectorXd::Zero(_nv_a * _no_a);
  Eigen::VectorXd nullVector_b = Eigen::VectorXd::Zero(_nv_b * _no_b);
  sigmaVector.head(_alpha) += this->getJ2GContribution(Zia_a.data(), nullptr, nullVector_a, true, 1);
  sigmaVector.tail(_beta) += this->getJ2GContribution(Zia_b.data(), nullptr, nullVector_b, true, -1);

  return sigmaVector;
}

template<>
Eigen::VectorXd CC2Controller<UNRESTRICTED>::getFDoubles(Eigen::Ref<Eigen::VectorXd> amp,
                                                         Eigen::Ref<Eigen::VectorXd> dip, double frequency) {
  // Dipole Integrals
  unsigned mo_a = _no_a + _nv_a;
  unsigned mo_b = _no_b + _nv_b;
  Eigen::Map<Eigen::MatrixXd> D_a(dip.data(), mo_a, mo_a);
  Eigen::Map<Eigen::MatrixXd> D_b(dip.data() + mo_a * mo_a, mo_b, mo_b);

  SPMatrix<UNRESTRICTED> dip_oo;
  SPMatrix<UNRESTRICTED> dip_vv;

  dip_oo.alpha = D_a.topLeftCorner(_no_a, _no_a);
  dip_vv.alpha = D_a.bottomRightCorner(_nv_a, _nv_a);

  dip_oo.beta = D_b.topLeftCorner(_no_b, _no_b);
  dip_vv.beta = D_b.bottomRightCorner(_nv_b, _nv_b);

  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);
  Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

  this->performTransformation(_Xia.alpha, amp.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, amp.tail(_beta), false, -1);
  _Zia.alpha.setZero();
  _Zia.beta.setZero();
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, frequency, 1) + this->getXiDoubles(i, j, frequency, dip_vv, 1);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
      Tia.segment(i * _nv_a, _nv_a).noalias() += rij * _gsLagrange.segment(j * _nv_a, _nv_a);
      if (i != j) {
        _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        Tia.segment(j * _nv_a, _nv_a).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_a, _nv_a);
      }
    }
  }
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, frequency, -1) + this->getXiDoubles(i, j, frequency, dip_vv, -1);
      _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
      Tia.segment(_alpha + i * _nv_b, _nv_b).noalias() += rij * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b);
      if (i != j) {
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
        Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(_alpha + i * _nv_b, _nv_b);
      }
    }
  }
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, frequency, 0) + this->getXiDoubles(i, j, frequency, dip_vv, 0);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
      Tia.segment(i * _nv_a, _nv_a).noalias() += rij * _gsLagrange.segment(_alpha + j * _nv_b, _nv_b);
      _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_a, _nv_a);
    }
  }

  this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
  if (_Bia != _Jia) {
    this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
  }
  this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
  this->reorderTensor(Tia.head(_alpha), _nv_a, _no_a);
  this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);

  this->reorderTensor(_Jia->beta, _nv_b, _no_b);
  if (_Bia != _Jia) {
    this->reorderTensor(_Bia->beta, _nv_b, _no_b);
  }
  this->reorderTensor(_Zia.beta, _nv_b, _no_b);
  this->reorderTensor(Tia.tail(_beta), _nv_b, _no_b);
  this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);

  for (size_t a = 0; a < _nv_a; ++a) {
    for (size_t b = a; b < _nv_a; ++b) {
      Eigen::MatrixXd rab = this->getXiDoublesV(a, b, frequency, dip_oo, 1);
      _Zia.alpha.middleRows(a * _no_a, _no_a).noalias() += rab * _Jia->alpha.middleRows(b * _no_a, _no_a);
      Tia.segment(a * _no_a, _no_a).noalias() += rab * _gsLagrange.segment(b * _no_a, _no_a);
      if (a != b) {
        _Zia.alpha.middleRows(b * _no_a, _no_a).noalias() += rab.transpose() * _Jia->alpha.middleRows(a * _no_a, _no_a);
        Tia.segment(b * _no_a, _no_a).noalias() += rab.transpose() * _gsLagrange.segment(a * _no_a, _no_a);
      }
    }
  }
  for (size_t a = 0; a < _nv_b; ++a) {
    for (size_t b = a; b < _nv_b; ++b) {
      Eigen::MatrixXd rab = this->getXiDoublesV(a, b, frequency, dip_oo, -1);
      _Zia.beta.middleRows(a * _no_b, _no_b).noalias() += rab * _Jia->beta.middleRows(b * _no_b, _no_b);
      Tia.segment(_alpha + a * _no_b, _no_b).noalias() += rab * _gsLagrange.segment(_alpha + b * _no_b, _no_b);
      if (a != b) {
        _Zia.beta.middleRows(b * _no_b, _no_b).noalias() += rab.transpose() * _Jia->beta.middleRows(a * _no_b, _no_b);
        Tia.segment(_alpha + b * _no_b, _no_b).noalias() += rab.transpose() * _gsLagrange.segment(_alpha + a * _no_b, _no_b);
      }
    }
  }
  for (size_t a = 0; a < _nv_a; ++a) {
    for (size_t b = 0; b < _nv_b; ++b) {
      Eigen::MatrixXd rab = this->getXiDoublesV(a, b, frequency, dip_oo, 0);
      _Zia.alpha.middleRows(a * _no_a, _no_a).noalias() += rab * _Jia->beta.middleRows(b * _no_b, _no_b);
      Tia.segment(a * _no_a, _no_a).noalias() += rab * _gsLagrange.segment(_alpha + b * _no_b, _no_b);
      _Zia.beta.middleRows(b * _no_b, _no_b).noalias() += rab.transpose() * _Jia->alpha.middleRows(a * _no_a, _no_a);
      Tia.segment(_alpha + b * _no_b, _no_b).noalias() += rab.transpose() * _gsLagrange.segment(a * _no_a, _no_a);
    }
  }

  this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
  if (_Bia != _Jia) {
    this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
  }
  this->reorderTensor(_Zia.alpha, _no_a, _nv_a);
  this->reorderTensor(Tia.head(_alpha), _no_a, _nv_a);
  this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);

  this->reorderTensor(_Jia->beta, _no_b, _nv_b);
  if (_Bia != _Jia) {
    this->reorderTensor(_Bia->beta, _no_b, _nv_b);
  }
  this->reorderTensor(_Zia.beta, _no_b, _nv_b);
  this->reorderTensor(Tia.tail(_beta), _no_b, _nv_b);
  this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  sigmaVector.head(_alpha) = _Jia->alpha * XQ;
  sigmaVector.tail(_beta) = _Jia->beta * XQ;

  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);

  // GH contribution
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    sigma_a.noalias() -= jia_a * tia_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    sigma_b.noalias() -= jia_b * tia_b.transpose() * jia_b;
  }

  return sigmaVector;
}

template<>
std::pair<std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>>
CC2Controller<UNRESTRICTED>::getFContractions(std::vector<Eigen::MatrixXd>& electricAmps,
                                              std::vector<Eigen::MatrixXd>& magneticAmps,
                                              std::vector<double> frequencies, Eigen::Ref<Eigen::MatrixXd> electricDips,
                                              Eigen::Ref<Eigen::MatrixXd> magneticDips, Options::GAUGE gauge) {
  unsigned nFreqs = frequencies.size();
  std::pair<std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>> dotProducts;

  dotProducts.first = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());
  dotProducts.second = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());

  double gaugeFactorForORD = (gauge == Options::GAUGE::VELOCITY) ? 1 : -1;

  auto dressIntegrals = [&](Eigen::Ref<Eigen::VectorXd> amp) {
    Timings::takeTime("CC2 -              Dress Ints");
    double dot = 0;
    Eigen::Map<Eigen::MatrixXd> t_a(amp.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
    Eigen::MatrixXd tg_a = t_a.transpose() * gsl_a;
    Eigen::MatrixXd gt_a = gsl_a * t_a.transpose();

    Eigen::Map<Eigen::MatrixXd> t_b(amp.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::MatrixXd tg_b = t_b.transpose() * gsl_b;
    Eigen::MatrixXd gt_b = gsl_b * t_b.transpose();

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::MatrixXd tmp_a = jia_a * tg_a + gt_a * jia_a;
      dot += zia_a.cwiseProduct(tmp_a).sum();

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
      Eigen::MatrixXd tmp_b = jia_b * tg_b + gt_b * jia_b;
      dot += zia_b.cwiseProduct(tmp_b).sum();
    }
    Timings::timeTaken("CC2 -              Dress Ints");
    return dot;
  };

  Timings::takeTime("CC2 -         F-Matrix Contr.");
  for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
    double frequency = frequencies[iFreq];
    unsigned p = iFreq;
    unsigned m = (frequencies[iFreq] != 0) ? iFreq + nFreqs : iFreq;
    double F = (frequencies[iFreq] == 0) ? 2.0 : 1.0;
    for (unsigned I = 0; I < 3; ++I) {
      Eigen::Ref<Eigen::VectorXd> epI = electricAmps[p].col(I);
      Eigen::Ref<Eigen::VectorXd> emI = electricAmps[m].col(I);

      Eigen::Ref<Eigen::VectorXd> mpI = magneticAmps[p].col(I);
      Eigen::Ref<Eigen::VectorXd> mmI = magneticAmps[m].col(I);

      /***********/
      /* Singles */
      /***********/
      Eigen::VectorXd eps = this->getFSingles(epI);
      // <r+|r->
      dotProducts.first[iFreq].row(I) += F * eps.transpose() * electricAmps[m];
      // <r+|m->
      dotProducts.second[iFreq].row(I) += gaugeFactorForORD * F * eps.transpose() * magneticAmps[m];

      if (frequency != 0) {
        Eigen::VectorXd ems = this->getFSingles(emI);
        // <r-|r+>
        dotProducts.first[iFreq].row(I) += ems.transpose() * electricAmps[p];
        // <r-|m+>
        dotProducts.second[iFreq].row(I) += ems.transpose() * magneticAmps[p];
      }

      /***********/
      /* Doubles */
      /***********/
      // Electric amplitudes plus frequency.
      auto epd = this->getFDoubles(epI, electricDips.col(I), frequency);
      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> emJ = electricAmps[m].col(J);
        double dot = epd.dot(emJ) - dressIntegrals(emJ);
        // <r+|r->
        dotProducts.first[iFreq](I, J) += F * dot;
        dotProducts.first[iFreq](J, I) += F * dot;
      }

      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> mmJ = magneticAmps[m].col(J);
        double dot = epd.dot(mmJ) - dressIntegrals(mmJ);
        // <r+|m->
        dotProducts.second[iFreq](I, J) += gaugeFactorForORD * F * dot;
      }

      // Magnetic amplitudes plus frequency.
      auto mpd = this->getFDoubles(mpI, magneticDips.col(I), frequency);
      for (unsigned J = 0; J < 3; ++J) {
        Eigen::Ref<Eigen::VectorXd> emJ = electricAmps[m].col(J);
        double dot = mpd.dot(emJ) - dressIntegrals(emJ);
        // <r-|m+>
        dotProducts.second[iFreq](J, I) += F * dot;
      }

      if (frequency != 0) {
        // Electric amplitudes minus frequency.
        auto emd = this->getFDoubles(emI, electricDips.col(I), -frequency);
        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> epJ = electricAmps[p].col(J);
          double dot = emd.dot(epJ) - dressIntegrals(epJ);
          // <r-|r+>
          dotProducts.first[iFreq](I, J) += dot;
          dotProducts.first[iFreq](J, I) += dot;
        }

        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> mpJ = magneticAmps[p].col(J);
          double dot = emd.dot(mpJ) - dressIntegrals(mpJ);
          // <r-|m+>
          dotProducts.second[iFreq](I, J) += dot;
        }

        // Magnetic amplitudes minus frequency.
        auto mmd = this->getFDoubles(mmI, magneticDips.col(I), -frequency);
        for (unsigned J = 0; J < 3; ++J) {
          Eigen::Ref<Eigen::VectorXd> epJ = electricAmps[p].col(J);
          double dot = mmd.dot(epJ) - dressIntegrals(epJ);
          // <r+|m->
          dotProducts.second[iFreq](J, I) += gaugeFactorForORD * dot;
        }
      }
    }
  } /* for iFreq */
  Timings::timeTaken("CC2 -         F-Matrix Contr.");

  return dotProducts;
} /* this->calcFContractions() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity