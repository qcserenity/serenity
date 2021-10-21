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
#include "postHF/LRSCF/Sigmavectors/RICC2/ADC2Sigmavector.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ADC2Sigmavector<SCFMode>::ADC2Sigmavector(std::shared_ptr<LRSCFController<SCFMode>> lrscf)
  : XWFController<SCFMode>(lrscf) {
}

template<>
void ADC2Sigmavector<RESTRICTED>::calculateE() {
  _Eij = Eigen::MatrixXd::Zero(_no, _no);
  _Eab = Eigen::MatrixXd::Zero(_nv, _nv);

  _Eij.diagonal() = _e.segment(0, _no);
  _Eab.diagonal() = _e.segment(_no, _nv);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.data() + Q * _nv * _no, _nv, _no);
    _Eij += jia.transpose() * yia;
    _Eab -= yia * jia.transpose();
  }

  // symmetrize
  _Eij = 0.5 * (_Eij + _Eij.transpose()).eval();
  _Eab = 0.5 * (_Eab + _Eab.transpose()).eval();

  // don't need Yia anymore for now
  _Yia.resize(0, 0);

  // gonna need these though
  _Xia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
  _Zia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
} /* this->calculateE() restricted */

template<>
void ADC2Sigmavector<UNRESTRICTED>::calculateE() {
  auto& Jia = *_Jia;
  for_spin(_Eij, _Eab, _nv, _no, _Yia, Jia, _e) {
    _Eij_spin = Eigen::MatrixXd::Zero(_no_spin, _no_spin);
    _Eab_spin = Eigen::MatrixXd::Zero(_nv_spin, _nv_spin);

    _Eij_spin.diagonal() = _e_spin.segment(0, _no_spin);
    _Eab_spin.diagonal() = _e_spin.segment(_no_spin, _nv_spin);

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(Jia_spin.data() + Q * _nv_spin * _no_spin, _nv_spin, _no_spin);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia_spin.data() + Q * _nv_spin * _no_spin, _nv_spin, _no_spin);
      _Eij_spin += jia.transpose() * yia;
      _Eab_spin -= yia * jia.transpose();
    }

    // symmetrize
    _Eij_spin = 0.5 * (_Eij_spin + _Eij_spin.transpose()).eval();
    _Eab_spin = 0.5 * (_Eab_spin + _Eab_spin.transpose()).eval();
  };

  // don't need Yia anymore now

  // gonna need these though
  _Xia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  _Zia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  for_spin(_Xia, _Yia, _Zia, _nv, _no) {
    _Yia_spin.resize(0, 0);
    _Xia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    _Zia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
  };
} /* this->calculateE() unrestricted */

template<>
Eigen::VectorXd ADC2Sigmavector<RESTRICTED>::getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia, guessVector, false);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  Eigen::VectorXd sigmaVector = 2 * _Jia->operator*(_Jia->transpose() * guessVector);
  Eigen::VectorXd Fia = sigmaVector;
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_no * _nv);

  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(Fia.data(), _nv, _no);

  sigma.noalias() += _Eab * guess - guess * _Eij;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _no * _nv, _nv, _no);
    fia.noalias() -= jia * guess.transpose() * jia;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Zia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j);
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
      Tia.segment(i * _nv, _nv).noalias() += tij * guessVector.segment(j * _nv, _nv);
      sigmaVector.segment(i * _nv, _nv).noalias() += 0.5 * tij * Fia.segment(j * _nv, _nv);
      _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
      if (i != j) {
        Tia.segment(j * _nv, _nv).noalias() += tij.transpose() * guessVector.segment(i * _nv, _nv);
        sigmaVector.segment(j * _nv, _nv).noalias() += 0.5 * tij.transpose() * Fia.segment(i * _nv, _nv);
        _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
      }
    }
  }
  Timings::timeTaken("Exc. State WF -  Amp-Contraction");

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  sigmaVector.noalias() += 0.5 * 2 * _Jia->operator*(_Jia->transpose() * Tia);
  Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij(_Jij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.data() + Q * _nv * _no, _nv, _no);
    sigma.noalias() -= 0.5 * jia * tia.transpose() * jia;
    sigma.noalias() -= zia * jij;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector.noalias() += this->getJ2GContribution(_Zia.data(), _Jij->data(), guessVector, false);

  return sigmaVector;
} /* this->getRightXWFSigma() restricted */

template<>
Eigen::VectorXd ADC2Sigmavector<UNRESTRICTED>::getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia.alpha, guessVector.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, guessVector.tail(_beta), false, -1);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
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
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _no_a * _nv_a, _nv_a, _no_a);
    fia_a.noalias() -= jia_a * guess_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _no_b * _nv_b, _nv_b, _no_b);
    fia_b.noalias() -= jia_b * guess_b.transpose() * jia_b;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Zia.alpha.setZero();
  _Zia.beta.setZero();
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

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
  sigmaVector.head(_alpha).noalias() += 0.5 * _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += 0.5 * _Jia->beta * XQ;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij_a(_Jij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    sigma_a.noalias() -= 0.5 * jia_a * tia_a.transpose() * jia_a;
    sigma_a.noalias() -= zia_a * jij_a;

    Eigen::Map<Eigen::MatrixXd> jij_b(_Jij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    sigma_b.noalias() -= 0.5 * jia_b * tia_b.transpose() * jia_b;
    sigma_b.noalias() -= zia_b * jij_b;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Zia.alpha.data(), _Jij->alpha.data(), guessVector.head(_alpha), false, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Zia.beta.data(), _Jij->beta.data(), guessVector.tail(_beta), false, -1);

  return sigmaVector;
} /* this->getRightXWFSigma() unrestricted */

template<>
void ADC2Sigmavector<RESTRICTED>::calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                                      std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("Exc. State WF - Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo = _no + _nv;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo * mo, nEigen);
  printf("  calculating %3i xi density matrices.\n\n", nEigen);

  for (size_t iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> D(densityMatrices[0].col(iEigen).data(), mo, mo);

    Eigen::Ref<Eigen::MatrixXd> Doo = D.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> Dov = D.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> Dvo = D.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> Dvv = D.bottomRightCorner(_nv, _nv);

    Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::Map<Eigen::MatrixXd> rev(eigenvector.data(), _nv, _no);

    this->performTransformation(_Xia, eigenvector, false);

    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudesA(i, j);
        Eigen::MatrixXd ttij = _soss * tij - _sss * tij.transpose();
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
        Dvv.noalias() += rij * tij.transpose();
        Dov.row(i).noalias() += ttij * rev.col(j);
        if (i != j) {
          Dvv.noalias() += rij.transpose() * tij;
          Dov.row(j).noalias() += ttij.transpose() * rev.col(i);
        }
      }
    }
    Dvo.noalias() += rev;

    // done with ij-wise amps, must reorder ints
    this->reorderTensor(*_Jia, _nv, _no, _nx);
    this->reorderTensor(_Xia, _nv, _no, _nx);

    for (size_t a = 0; a < _nv; ++a) {
      for (size_t b = a; b < _nv; ++b) {
        Eigen::MatrixXd tab = this->getAmplitudesAV(a, b);
        Eigen::MatrixXd rab = this->getRightAmplitudesV(a, b, eigenvalue);
        Doo.noalias() -= tab * rab.transpose();
        if (a != b) {
          Doo.noalias() -= tab.transpose() * rab;
        }
      }
    }

    // done with ab-wise amps, must order ints back
    this->reorderTensor(*_Jia, _no, _nv, _nx);
  } /* iEigen */

  // copy right into left density matrices
  densityMatrices[1] = densityMatrices[0];

  Timings::timeTaken("Exc. State WF - Density Matrices");
} /* this->calculateDensityMatrices() restricted */

template<>
void ADC2Sigmavector<UNRESTRICTED>::calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                                        std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("Exc. State WF - Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned alphaMO = _no_a + _nv_a;
  unsigned betaMO = _no_b + _nv_b;
  densityMatrices[0] = Eigen::MatrixXd::Zero(alphaMO * alphaMO + betaMO * betaMO, nEigen);
  printf("  calculating %3i xi density matrices.\n\n", nEigen);

  for (size_t iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[0].col(iEigen).data(), alphaMO, alphaMO);
    Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[0].col(iEigen).data() + alphaMO * alphaMO, betaMO, betaMO);

    Eigen::Ref<Eigen::MatrixXd> Doo_a = D_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> Dov_a = D_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> Dvo_a = D_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> Dvv_a = D_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> Doo_b = D_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> Dov_b = D_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> Dvo_b = D_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> Dvv_b = D_b.bottomRightCorner(_nv_b, _nv_b);

    Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::Map<Eigen::MatrixXd> rev_a(eigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> rev_b(eigenvector.data() + _alpha, _nv_b, _no_b);

    this->performTransformation(_Xia.alpha, eigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, eigenvector.tail(_beta), false, -1);

    // alpha
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudesA(i, j, 1);
        Eigen::MatrixXd ttij = _sss * tij - _sss * tij.transpose();
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 1);
        Dvv_a.noalias() += rij * tij.transpose();
        Dov_a.row(i).noalias() += ttij * rev_a.col(j);
        if (i != j) {
          Dvv_a.noalias() += rij.transpose() * tij;
          Dov_a.row(j).noalias() += ttij.transpose() * rev_a.col(i);
        }
      }
    }

    // beta
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudesA(i, j, -1);
        Eigen::MatrixXd ttij = _sss * tij - _sss * tij.transpose();
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
        Dvv_b.noalias() += rij * tij.transpose();
        Dov_b.row(i).noalias() += ttij * rev_b.col(j);
        if (i != j) {
          Dvv_b.noalias() += rij.transpose() * tij;
          Dov_b.row(j).noalias() += ttij.transpose() * rev_b.col(i);
        }
      }
    }

    // mixed
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getAmplitudesA(i, j, 0);
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);
        Dvv_a.noalias() += rij * tij.transpose();
        Dov_a.row(i).noalias() += _oss * tij * rev_b.col(j);
        Dvv_b.noalias() += rij.transpose() * tij;
        Dov_b.row(j).noalias() += _oss * tij.transpose() * rev_a.col(i);
      }
    }
    Dvo_a += rev_a;
    Dvo_b += rev_b;

    // done with ij-wise amps, must reorder ints
    this->reorderTensor(_Jia->alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Xia.alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Jia->beta, _nv_b, _no_b, _nx);
    this->reorderTensor(_Xia.beta, _nv_b, _no_b, _nx);

    // alpha
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = a; b < _nv_a; ++b) {
        Eigen::MatrixXd tab = this->getAmplitudesAV(a, b, 1);
        Eigen::MatrixXd rab = this->getRightAmplitudesV(a, b, eigenvalue, 1);
        Doo_a.noalias() -= tab * rab.transpose();
        if (a != b) {
          Doo_a.noalias() -= tab.transpose() * rab;
        }
      }
    }

    // beta
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd tab = this->getAmplitudesAV(a, b, -1);
        Eigen::MatrixXd rab = this->getRightAmplitudesV(a, b, eigenvalue, -1);
        Doo_b.noalias() -= tab * rab.transpose();
        if (a != b) {
          Doo_b.noalias() -= tab.transpose() * rab;
        }
      }
    }

    // mixed
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd tab = this->getAmplitudesAV(a, b, 0);
        Eigen::MatrixXd rab = this->getRightAmplitudesV(a, b, eigenvalue, 0);
        Doo_a.noalias() -= tab * rab.transpose();
        Doo_b.noalias() -= tab.transpose() * rab;
      }
    }

    // done with ab-wise amps, must order ints back
    this->reorderTensor(_Jia->alpha, _no_a, _nv_a, _nx);
    this->reorderTensor(_Jia->beta, _no_b, _nv_b, _nx);
  } /* iEigen */

  // copy right into left density matrices
  densityMatrices[1] = densityMatrices[0];

  Timings::timeTaken("Exc. State WF - Density Matrices");
} /* this->calculateDensityMatrices() unrestricted */

template class ADC2Sigmavector<Options::SCF_MODES::RESTRICTED>;
template class ADC2Sigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity