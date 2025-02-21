/**
 * @file CC2DensityTransition.cpp
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
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "io/FormattedOutputStream.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::calculateTransitionDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                             Eigen::VectorXd eigenvalues,
                                                             std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo = _no + _nv;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo * mo, nEigen);
  densityMatrices[1] = Eigen::MatrixXd::Zero(mo * mo, nEigen);

  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, 2 * nEigen);

  Eigen::MatrixXd CC2virt = 0.5 * _gsDensity.bottomRightCorner(_nv, _nv);
  Eigen::MatrixXd CC2occ = 0.5 * _gsDensity.topLeftCorner(_no, _no);
  CC2occ.diagonal() -= Eigen::VectorXd::Constant(_no, 1);

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    OutputControl::v.printf("  Calculating transition density matrices for state %3i of %3i.\n", iEigen + 1, nEigen);

    Eigen::Map<Eigen::MatrixXd> rightD(densityMatrices[0].col(iEigen).data(), mo, mo);
    Eigen::Map<Eigen::MatrixXd> leftD(densityMatrices[1].col(iEigen).data(), mo, mo);

    Eigen::Ref<Eigen::MatrixXd> rightDoo = rightD.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> rightDov = rightD.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> rightDvo = rightD.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> rightDvv = rightD.bottomRightCorner(_nv, _nv);

    Eigen::Ref<Eigen::MatrixXd> leftDoo = leftD.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> leftDov = leftD.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> leftDvo = leftD.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> leftDvv = leftD.bottomRightCorner(_nv, _nv);

    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[3].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);

    Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> rev(rightEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> lev(leftEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> exl(excLagrange.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fai(_Fai.data(), _nv, _no);

    _Fai = 2 * _Jia->operator*(_Jia->transpose() * rightEigenvector);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      fai.noalias() -= jia * rev.transpose() * jia;
    }

    _leftSingles = leftEigenvector;
    _exlSingles1 = excLagrange;
    _exlSingles2 = _gsLagrange;

    // Easy contributions first
    rightDov.noalias() += rev.transpose();
    rightDvo.noalias() += exl;
    leftDvo.noalias() += lev;
    rightDov.noalias() -= rev.transpose() * CC2virt;
    rightDov.noalias() += CC2occ * rev.transpose();
    rightDvv.noalias() += gsl * rev.transpose();
    rightDoo.noalias() -= rev.transpose() * gsl;

    this->performTransformation(_Xia, rightEigenvector, false);
    this->performTransformation(_Zia, _gsLagrange, true);
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
        Eigen::MatrixXd rrij = _soss * rij - _sss * rij.transpose();
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        rightDvv.noalias() += tij * rij.transpose();
        rightDov.row(i).noalias() += rrij * gsl.col(j);
        if (i != j) {
          rightDvv.noalias() += tij.transpose() * rij;
          rightDov.row(j).noalias() += rrij.transpose() * gsl.col(i);
        }
      }
    }

    this->performTransformation(_Yia, excLagrange, true);
    this->performTransformation(_Zia, leftEigenvector, true);
    Eigen::MatrixXd rg = -1 * rev.transpose() * gsl;
    Eigen::MatrixXd gr = -1 * gsl * rev.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> wia(_Wia.col(Q).data(), _nv, _no);
      wia = jia * rg + gr * jia;
    }

    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j);
        Eigen::MatrixXd iij = _soss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue);
        leftDvv.noalias() += lij * ij.transpose();
        leftDov.row(i).noalias() += iij * lev.col(j);
        rightDvv.noalias() += mij * ij.transpose();
        rightDov.row(i).noalias() += iij * exl.col(j);
        if (i != j) {
          leftDvv.noalias() += lij.transpose() * ij;
          leftDov.row(j).noalias() += iij.transpose() * lev.col(i);
          rightDvv.noalias() += mij.transpose() * ij;
          rightDov.row(j).noalias() += iij.transpose() * exl.col(i);
        }
      }
    }

    // Done with ij-wise amplitudes, must reorder integrals.
    this->reorderTensor(*_Jia, _nv, _no);
    if (_Jia != _Bia) {
      this->reorderTensor(*_Bia, _nv, _no);
    }
    this->reorderTensor(_Fia, _nv, _no);
    this->reorderTensor(_gsLagrange, _nv, _no);

    this->reorderTensor(_Wia, _nv, _no);
    this->reorderTensor(_Xia, _nv, _no);
    this->reorderTensor(_Yia, _nv, _no);
    this->reorderTensor(_Zia, _nv, _no);

    this->reorderTensor(_Fai, _nv, _no);

    this->reorderTensor(_leftSingles, _nv, _no);
    this->reorderTensor(_exlSingles1, _nv, _no);
    this->reorderTensor(_exlSingles2, _nv, _no);

    for (size_t a = 0; a < _nv; ++a) {
      for (size_t b = a; b < _nv; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue);
        leftDoo.noalias() -= ab * lab.transpose();
        rightDoo.noalias() -= ab * mab.transpose();
        if (a != b) {
          leftDoo.noalias() -= ab.transpose() * lab;
          rightDoo.noalias() -= ab.transpose() * mab;
        }
      }
    }
    this->reorderTensor(_gsLagrange, _no, _nv);
    this->performTransformation(_Zia, _gsLagrange, true);
    this->reorderTensor(_gsLagrange, _nv, _no);
    this->reorderTensor(_Zia, _nv, _no);
    for (size_t a = 0; a < _nv; ++a) {
      for (size_t b = a; b < _nv; ++b) {
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue);
        rightDoo.noalias() -= rab * tab.transpose();
        if (a != b) {
          rightDoo.noalias() -= rab.transpose() * tab;
        }
      }
    }

    // Done with ab-wise amplitudes, must order integrals back.
    this->reorderTensor(*_Jia, _no, _nv);
    if (_Jia != _Bia) {
      this->reorderTensor(*_Bia, _no, _nv);
    }
    this->reorderTensor(_Fia, _no, _nv);
    this->reorderTensor(_gsLagrange, _no, _nv);
  } /* iEigen */

  Timings::timeTaken("CC2 -        Density Matrices");
}

template<>
void CC2Controller<UNRESTRICTED>::calculateTransitionDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                               Eigen::VectorXd eigenvalues,
                                                               std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo_a = _no_a + _nv_a;
  unsigned mo_b = _no_b + _nv_b;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo_a * mo_a + mo_b * mo_b, nEigen);
  densityMatrices[1] = Eigen::MatrixXd::Zero(mo_a * mo_a + mo_b * mo_b, nEigen);

  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, 2 * nEigen);

  Eigen::MatrixXd CC2virt_a = _gsDensity.alpha.bottomRightCorner(_nv_a, _nv_a);
  Eigen::MatrixXd CC2occ_a = _gsDensity.alpha.topLeftCorner(_no_a, _no_a);
  CC2occ_a.diagonal() -= Eigen::VectorXd::Constant(_no_a, 1);

  Eigen::MatrixXd CC2virt_b = _gsDensity.beta.bottomRightCorner(_nv_b, _nv_b);
  Eigen::MatrixXd CC2occ_b = _gsDensity.beta.topLeftCorner(_no_b, _no_b);
  CC2occ_b.diagonal() -= Eigen::VectorXd::Constant(_no_b, 1);

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> rightD_a(densityMatrices[0].col(iEigen).data(), mo_a, mo_a);
    Eigen::Map<Eigen::MatrixXd> leftD_a(densityMatrices[1].col(iEigen).data(), mo_a, mo_a);

    Eigen::Map<Eigen::MatrixXd> rightD_b(densityMatrices[0].col(iEigen).data() + mo_a * mo_a, mo_b, mo_b);
    Eigen::Map<Eigen::MatrixXd> leftD_b(densityMatrices[1].col(iEigen).data() + mo_a * mo_a, mo_b, mo_b);

    Eigen::Ref<Eigen::MatrixXd> rightDoo_a = rightD_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> rightDov_a = rightD_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> rightDvo_a = rightD_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> rightDvv_a = rightD_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> leftDoo_a = leftD_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> leftDov_a = leftD_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> leftDvo_a = leftD_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> leftDvv_a = leftD_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> rightDoo_b = rightD_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> rightDov_b = rightD_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> rightDvo_b = rightD_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> rightDvv_b = rightD_b.bottomRightCorner(_nv_b, _nv_b);

    Eigen::Ref<Eigen::MatrixXd> leftDoo_b = leftD_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> leftDov_b = leftD_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> leftDvo_b = leftD_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> leftDvv_b = leftD_b.bottomRightCorner(_nv_b, _nv_b);

    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[3].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);

    Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> rev_a(rightEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> lev_a(leftEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> exl_a(excLagrange.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fai_a(_Fai.data(), _nv_a, _no_a);

    Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> rev_b(rightEigenvector.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> lev_b(leftEigenvector.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> exl_b(excLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fai_b(_Fai.data() + _alpha, _nv_b, _no_b);

    Eigen::VectorXd XQ =
        _Jia->alpha.transpose() * rightEigenvector.head(_alpha) + _Jia->beta.transpose() * rightEigenvector.tail(_beta);
    _Fai.head(_alpha) = _Jia->alpha * XQ;
    _Fai.tail(_beta) = _Jia->beta * XQ;
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      fai_a.noalias() -= jia_a * rev_a.transpose() * jia_a;
      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      fai_b.noalias() -= jia_b * rev_b.transpose() * jia_b;
    }

    _leftSingles = leftEigenvector;
    _exlSingles1 = excLagrange;
    _exlSingles2 = _gsLagrange;

    // Easy contributions first
    rightDov_a.noalias() += rev_a.transpose();
    rightDvo_a.noalias() += exl_a;
    leftDvo_a.noalias() += lev_a;
    rightDov_a.noalias() -= rev_a.transpose() * CC2virt_a;
    rightDov_a.noalias() += CC2occ_a * rev_a.transpose();
    rightDvv_a.noalias() += gsl_a * rev_a.transpose();
    rightDoo_a.noalias() -= rev_a.transpose() * gsl_a;

    rightDov_b.noalias() += rev_b.transpose();
    rightDvo_b.noalias() += exl_b;
    leftDvo_b.noalias() += lev_b;
    rightDov_b.noalias() -= rev_b.transpose() * CC2virt_b;
    rightDov_b.noalias() += CC2occ_b * rev_b.transpose();
    rightDvv_b.noalias() += gsl_b * rev_b.transpose();
    rightDoo_b.noalias() -= rev_b.transpose() * gsl_b;

    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        rightDvv_a.noalias() += tij * rij.transpose();
        rightDov_a.row(i).noalias() += rrij * gsl_a.col(j);
        if (i != j) {
          rightDvv_a.noalias() += tij.transpose() * rij;
          rightDov_a.row(j).noalias() += rrij.transpose() * gsl_a.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        rightDvv_b.noalias() += tij * rij.transpose();
        rightDov_b.row(i).noalias() += rrij * gsl_b.col(j);
        if (i != j) {
          rightDvv_b.noalias() += tij.transpose() * rij;
          rightDov_b.row(j).noalias() += rrij.transpose() * gsl_b.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        rightDvv_a.noalias() += tij * rij.transpose();
        rightDov_a.row(i).noalias() += _oss * rij * gsl_b.col(j);
        rightDvv_b.noalias() += tij.transpose() * rij;
        rightDov_b.row(j).noalias() += _oss * rij.transpose() * gsl_a.col(i);
      }
    }

    this->performTransformation(_Yia.alpha, excLagrange.head(_alpha), true, 1);
    this->performTransformation(_Yia.beta, excLagrange.tail(_beta), true, -1);
    this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);
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

    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, 1);
        leftDvv_a.noalias() += lij * ij.transpose();
        leftDov_a.row(i).noalias() += iij * lev_a.col(j);
        rightDvv_a.noalias() += mij * ij.transpose();
        rightDov_a.row(i).noalias() += iij * exl_a.col(j);
        if (i != j) {
          leftDvv_a.noalias() += lij.transpose() * ij;
          leftDov_a.row(j).noalias() += iij.transpose() * lev_a.col(i);
          rightDvv_a.noalias() += mij.transpose() * ij;
          rightDov_a.row(j).noalias() += iij.transpose() * exl_a.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, -1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, -1);
        leftDvv_b.noalias() += lij * ij.transpose();
        leftDov_b.row(i).noalias() += iij * lev_b.col(j);
        rightDvv_b.noalias() += mij * ij.transpose();
        rightDov_b.row(i).noalias() += iij * exl_b.col(j);
        if (i != j) {
          leftDvv_b.noalias() += lij.transpose() * ij;
          leftDov_b.row(j).noalias() += iij.transpose() * lev_b.col(i);
          rightDvv_b.noalias() += mij.transpose() * ij;
          rightDov_b.row(j).noalias() += iij.transpose() * exl_b.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 0);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, 0);
        leftDvv_a.noalias() += lij * ij.transpose();
        leftDov_a.row(i).noalias() += _oss * ij * lev_b.col(j);
        rightDvv_a.noalias() += mij * ij.transpose();
        rightDov_a.row(i).noalias() += _oss * ij * exl_b.col(j);
        leftDvv_b.noalias() += lij.transpose() * ij;
        leftDov_b.row(j).noalias() += _oss * ij.transpose() * lev_a.col(i);
        rightDvv_b.noalias() += mij.transpose() * ij;
        rightDov_b.row(j).noalias() += _oss * ij.transpose() * exl_a.col(i);
      }
    }

    // Done with ij-wise amplitudes, must reorder integrals.
    this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
    this->reorderTensor(_Jia->beta, _nv_b, _no_b);
    if (_Jia != _Bia) {
      this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
      this->reorderTensor(_Bia->beta, _nv_b, _no_b);
    }
    this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
    this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);

    this->reorderTensor(_Wia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Wia.beta, _nv_b, _no_b);
    this->reorderTensor(_Xia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Xia.beta, _nv_b, _no_b);
    this->reorderTensor(_Yia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Yia.beta, _nv_b, _no_b);
    this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Zia.beta, _nv_b, _no_b);

    this->reorderTensor(_Fai.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_Fai.tail(_beta), _nv_b, _no_b);

    this->reorderTensor(_leftSingles.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_leftSingles.tail(_beta), _nv_b, _no_b);
    this->reorderTensor(_exlSingles1.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_exlSingles1.tail(_beta), _nv_b, _no_b);
    this->reorderTensor(_exlSingles2.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_exlSingles2.tail(_beta), _nv_b, _no_b);

    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = a; b < _nv_a; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 1);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, 1);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, 1);
        leftDoo_a.noalias() -= ab * lab.transpose();
        rightDoo_a.noalias() -= ab * mab.transpose();
        if (a != b) {
          leftDoo_a.noalias() -= ab.transpose() * lab;
          rightDoo_a.noalias() -= ab.transpose() * mab;
        }
      }
    }
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, -1);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, -1);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, -1);
        leftDoo_b.noalias() -= ab * lab.transpose();
        rightDoo_b.noalias() -= ab * mab.transpose();
        if (a != b) {
          leftDoo_b.noalias() -= ab.transpose() * lab;
          rightDoo_b.noalias() -= ab.transpose() * mab;
        }
      }
    }
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 0);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, 0);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, 0);
        leftDoo_a.noalias() -= ab * lab.transpose();
        leftDoo_b.noalias() -= ab.transpose() * lab;
        rightDoo_a.noalias() -= ab * mab.transpose();
        rightDoo_b.noalias() -= ab.transpose() * mab;
      }
    }
    this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);
    this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Zia.beta, _nv_b, _no_b);
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = a; b < _nv_a; ++b) {
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 1);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 1);
        rightDoo_a.noalias() -= rab * tab.transpose();
        if (a != b) {
          rightDoo_a.noalias() -= rab.transpose() * tab;
        }
      }
    }
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, -1);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, -1);
        rightDoo_b.noalias() -= rab * tab.transpose();
        if (a != b) {
          rightDoo_b.noalias() -= rab.transpose() * tab;
        }
      }
    }
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 0);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 0);
        rightDoo_a.noalias() -= rab * tab.transpose();
        rightDoo_b.noalias() -= rab.transpose() * tab;
      }
    }

    // Done with ab-wise amplitudes, must order integrals back.
    this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
    this->reorderTensor(_Jia->beta, _no_b, _nv_b);
    if (_Jia != _Bia) {
      this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
      this->reorderTensor(_Bia->beta, _no_b, _nv_b);
    }
    this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
    this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
    this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
  } /* iEigen */

  Timings::timeTaken("CC2 -        Density Matrices");
}

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity