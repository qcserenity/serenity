/**
 * @file CC2MultiplierTransition.cpp
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
#include <iomanip>

namespace Serenity {

template<>
std::vector<Eigen::MatrixXd>
CC2Controller<RESTRICTED>::calculateTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                       Eigen::VectorXd eigenvalues) {
  unsigned nEigen = eigenvalues.size();
  std::vector<Eigen::MatrixXd> rightHandSides(nEigen, Eigen::MatrixXd::Zero(_no * _nv, 1));

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // get eigenpair data
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides[iEigen].col(0);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nv * _no);
    Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nv * _no);

    Eigen::Map<Eigen::MatrixXd> rev(rightEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> rhs(rightHandSide.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fai(_Fai.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

    Eigen::MatrixXd Eij = Eigen::MatrixXd::Zero(_no, _no);
    Eigen::MatrixXd Eab = Eigen::MatrixXd::Zero(_nv, _nv);

    this->performTransformation(_Xia, rightEigenvector, false);

    Eij += fia.transpose() * rev;
    Eab -= rev * fia.transpose();

    _Zia.setZero();
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
        _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
        Tia.segment(i * _nv, _nv).noalias() += rij * _gsLagrange.segment(j * _nv, _nv);
        if (i != j) {
          _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
          Tia.segment(j * _nv, _nv).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv, _nv);
        }
      }
    }

    Eigen::VectorXd XQ = _Jia->transpose() * rightEigenvector;
    _Fai = 2 * _Jia->operator*(XQ);

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
      // ^Fij and ^Fab
      Eij.noalias() += 2 * bij * XQ(Q);

      Eij.noalias() -= jia.transpose() * rev * bij;

      Eij.noalias() += jia.transpose() * zia;
      Eab.noalias() -= zia * jia.transpose();

      fai.noalias() -= jia * rev.transpose() * jia;
    }

    Eab += this->transformedABFock(rightEigenvector, XQ, _Yia);

    // 0 contribution
    rhs += Eab.transpose() * gsl - gsl * Eij.transpose();

    // I contribution
    Tia += rightEigenvector;
    rightHandSide += 2 * _Jia->operator*(_Jia->transpose() * Tia);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      rhs.noalias() -= jia * tia.transpose() * jia;
    }

    // GH contribution
    Eigen::MatrixXd rg = -1 * rev.transpose() * gsl;
    Eigen::MatrixXd gr = -1 * gsl * rev.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> wia(_Wia.col(Q).data(), _nv, _no);
      wia = jia * rg + gr * jia;
    }

    _Yia.setZero();
    _fockSingles = _gsLagrange;
    this->performTransformation(_Zia, _gsLagrange, true);
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue);
        _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Xia.middleRows(j * _nv, _nv);
        _Yia.middleRows(i * _nv, _nv).noalias() += fij * _Bia->middleRows(j * _nv, _nv);
        if (i != j) {
          _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Xia.middleRows(i * _nv, _nv);
          _Yia.middleRows(j * _nv, _nv).noalias() += fij.transpose() * _Bia->middleRows(i * _nv, _nv);
        }
      }
    }

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
      rhs -= yia * bij.transpose();
    }
    rightHandSide += this->getJ2GContribution(_Yia.data(), nullptr, nullVector, true);

    _Yia.setZero();
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Bia->middleRows(j * _nv, _nv);
        if (i != j) {
          _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Bia->middleRows(i * _nv, _nv);
        }
      }
    }

    rightHandSide += 2 * _Jia->operator*(_Xia.transpose() * _gsLagrange);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab = -1 * rev * jia.transpose();
      Eigen::MatrixXd _bij = jia.transpose() * rev;

      rhs.noalias() += _bab.transpose() * yia;
      rhs.noalias() -= yia * _bij.transpose();
      rhs.noalias() -= _bab.transpose() * gsl * bij.transpose();

      yia.noalias() = gsl * _bij.transpose();
    }
    rightHandSide -= this->getJ2GContribution(_Yia.data(), nullptr, nullVector, true);

    rightHandSide *= -1;
  }

  return rightHandSides;
} /* this->calculateTransitionMomentLagrangeMultiplier() restricted */

template<>
std::vector<Eigen::MatrixXd>
CC2Controller<UNRESTRICTED>::calculateTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                         Eigen::VectorXd eigenvalues) {
  unsigned nEigen = eigenvalues.size();
  std::vector<Eigen::MatrixXd> rightHandSides(nEigen, Eigen::MatrixXd::Zero(_nDim, 1));

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // get eigenpair data
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides[iEigen].col(0);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);
    Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nDim);
    _Fai = Eigen::VectorXd::Zero(_nDim);

    Eigen::Map<Eigen::MatrixXd> rev_a(rightEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> rhs_a(rightHandSide.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fai_a(_Fai.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);

    Eigen::Map<Eigen::MatrixXd> rev_b(rightEigenvector.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> rhs_b(rightHandSide.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fai_b(_Fai.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

    Eigen::MatrixXd Eij_a = Eigen::MatrixXd::Zero(_no_a, _no_a);
    Eigen::MatrixXd Eab_a = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
    Eigen::MatrixXd Eij_b = Eigen::MatrixXd::Zero(_no_b, _no_b);
    Eigen::MatrixXd Eab_b = Eigen::MatrixXd::Zero(_nv_b, _nv_b);

    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);

    Eij_a += fia_a.transpose() * rev_a;
    Eab_a -= rev_a * fia_a.transpose();
    Eij_b += fia_b.transpose() * rev_b;
    Eab_b -= rev_b * fia_b.transpose();

    _Zia.alpha.setZero();
    _Zia.beta.setZero();
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 1);
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
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
        _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        Tia.segment(i * _nv_b + _alpha, _nv_b).noalias() += rij * _gsLagrange.segment(j * _nv_b + _alpha, _nv_b);
        if (i != j) {
          _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
          Tia.segment(j * _nv_b + _alpha, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_b + _alpha, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        Tia.segment(i * _nv_a, _nv_a).noalias() += rij * _gsLagrange.segment(j * _nv_b + _alpha, _nv_b);
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        Tia.segment(j * _nv_b + _alpha, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_a, _nv_a);
      }
    }

    Eigen::VectorXd XQ =
        _Jia->alpha.transpose() * rightEigenvector.head(_alpha) + _Jia->beta.transpose() * rightEigenvector.tail(_beta);
    _Fai.head(_alpha) = _Jia->alpha * XQ;
    _Fai.tail(_beta) = _Jia->beta * XQ;

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
      // ^Fij and ^Fab
      Eij_a.noalias() += bij_a * XQ(Q);

      Eij_a.noalias() -= jia_a.transpose() * rev_a * bij_a;

      Eij_a.noalias() += jia_a.transpose() * zia_a;
      Eab_a.noalias() -= zia_a * jia_a.transpose();

      fai_a.noalias() -= jia_a * rev_a.transpose() * jia_a;

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
      // ^Fij and ^Fab
      Eij_b.noalias() += bij_b * XQ(Q);

      Eij_b.noalias() -= jia_b.transpose() * rev_b * bij_b;

      Eij_b.noalias() += jia_b.transpose() * zia_b;
      Eab_b.noalias() -= zia_b * jia_b.transpose();

      fai_b.noalias() -= jia_b * rev_b.transpose() * jia_b;
    }

    Eab_a += this->transformedABFock(rightEigenvector.head(_alpha), XQ, _Yia.alpha, 1);
    Eab_b += this->transformedABFock(rightEigenvector.tail(_beta), XQ, _Yia.beta, -1);

    // 0 contribution
    rhs_a += Eab_a.transpose() * gsl_a - gsl_a * Eij_a.transpose();
    rhs_b += Eab_b.transpose() * gsl_b - gsl_b * Eij_b.transpose();

    // I contribution
    Tia += rightEigenvector;
    XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
    rightHandSide.head(_alpha) += _Jia->alpha * XQ;
    rightHandSide.tail(_beta) += _Jia->beta * XQ;
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      rhs_a -= jia_a * tia_a.transpose() * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      rhs_b -= jia_b * tia_b.transpose() * jia_b;
    }

    // GH contribution
    Eigen::MatrixXd rg_a = -1 * rev_a.transpose() * gsl_a;
    Eigen::MatrixXd gr_a = -1 * gsl_a * rev_a.transpose();
    Eigen::MatrixXd rg_b = -1 * rev_b.transpose() * gsl_b;
    Eigen::MatrixXd gr_b = -1 * gsl_b * rev_b.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> wia_a(_Wia.alpha.col(Q).data(), _nv_a, _no_a);
      wia_a = jia_a * rg_a + gr_a * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> wia_b(_Wia.beta.col(Q).data(), _nv_b, _no_b);
      wia_b = jia_b * rg_b + gr_b * jia_b;
    }

    _Yia.alpha.setZero();
    _Yia.beta.setZero();
    _fockSingles = _gsLagrange;
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, 1);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.alpha.middleRows(j * _nv_a, _nv_a);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += fij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += fij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, -1);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += fij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.beta.middleRows(i * _nv_b, _nv_b);
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += fij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, 0);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += fij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += fij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);
      rhs_a -= yia_a * bij_a.transpose();

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);
      rhs_b -= yia_b * bij_b.transpose();
    }
    rightHandSide.head(_alpha) += this->getJ2GContribution(_Yia.alpha.data(), nullptr, nullVector, true, 1);
    rightHandSide.tail(_beta) += this->getJ2GContribution(_Yia.beta.data(), nullptr, nullVector, true, -1);

    _Yia.alpha.setZero();
    _Yia.beta.setZero();
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }

    XQ = _Xia.alpha.transpose() * _gsLagrange.head(_alpha) + _Xia.beta.transpose() * _gsLagrange.tail(_beta);
    rightHandSide.head(_alpha) += _Jia->alpha * XQ;
    rightHandSide.tail(_beta) += _Jia->beta * XQ;
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab_a = -1 * rev_a * jia_a.transpose();
      Eigen::MatrixXd _bij_a = jia_a.transpose() * rev_a;

      rhs_a.noalias() += _bab_a.transpose() * yia_a;
      rhs_a.noalias() -= yia_a * _bij_a.transpose();
      rhs_a.noalias() -= _bab_a.transpose() * gsl_a * bij_a.transpose();

      yia_a.noalias() = gsl_a * _bij_a.transpose();

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab_b = -1 * rev_b * jia_b.transpose();
      Eigen::MatrixXd _bij_b = jia_b.transpose() * rev_b;

      rhs_b.noalias() += _bab_b.transpose() * yia_b;
      rhs_b.noalias() -= yia_b * _bij_b.transpose();
      rhs_b.noalias() -= _bab_b.transpose() * gsl_b * bij_b.transpose();

      yia_b.noalias() = gsl_b * _bij_b.transpose();
    }
    rightHandSide.head(_alpha) -= this->getJ2GContribution(_Yia.alpha.data(), nullptr, nullVector, true, 1);
    rightHandSide.tail(_beta) -= this->getJ2GContribution(_Yia.beta.data(), nullptr, nullVector, true, -1);

    rightHandSide *= -1;
  }

  return rightHandSides;
} /* this->calculateTransitionMomentLagrangeMultiplier() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity