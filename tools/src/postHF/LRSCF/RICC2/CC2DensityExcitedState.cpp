/**
 * @file CC2DensityExcitedState.cpp
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
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::calculateExcitedStateDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                               Eigen::VectorXd eigenvalues,
                                                               std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");
  this->calculateGroundStateDensity();

  unsigned nEigen = _settings.ccexdens ? eigenvalues.size() : 0;
  unsigned mo = _no + _nv;
  densityMatrices[2] = Eigen::MatrixXd::Zero(mo * mo, 1 + nEigen);

  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, nEigen);

  // First column: only ground-state density.
  densityMatrices[2].colwise() += Eigen::Map<Eigen::VectorXd>(_gsDensity.data(), _gsDensity.size());

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> D(densityMatrices[2].col(1 + iEigen).data(), mo, mo);

    Eigen::MatrixXd virt = Eigen::MatrixXd::Zero(_nv, _nv);
    Eigen::MatrixXd occ = Eigen::MatrixXd::Zero(_no, _no);

    Eigen::Ref<Eigen::MatrixXd> Doo = D.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> Dov = D.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> Dvo = D.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> Dvv = D.bottomRightCorner(_nv, _nv);

    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[2].col(1 + iEigen);
    double eigenvalue = eigenvalues(iEigen);

    Eigen::Map<Eigen::MatrixXd> rev(rightEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> lev(leftEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> exl(excLagrange.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fai(_Fai.data(), _nv, _no);

    _Fai = 2 * _Jia->operator*(_Jia->transpose() * rightEigenvector);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      fai.noalias() -= jia * rev.transpose() * jia;
    }

    this->performTransformation(_Xia, rightEigenvector, false);
    this->performTransformation(_Yia, excLagrange, true);
    this->performTransformation(_Zia, leftEigenvector, true);

    Eigen::MatrixXd rg = -1 * rev.transpose() * lev;
    Eigen::MatrixXd gr = -1 * lev * rev.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
      Eigen::Map<Eigen::MatrixXd> wia(_Wia.col(Q).data(), _nv, _no);
      wia = jia * rg + gr * jia;
    }

    _leftSingles = leftEigenvector;
    _exlSingles1 = excLagrange;
    _exlSingles2 = leftEigenvector;

    Doo.noalias() -= rev.transpose() * lev;
    Dvv.noalias() += lev * rev.transpose();
    Dvo.noalias() += exl;

    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j);
        Eigen::MatrixXd iij = _soss * ij - _sss * ij.transpose();

        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
        Eigen::MatrixXd rrij = _soss * rij - _sss * rij.transpose();
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, 0);

        Dvv.noalias() += lij * rij.transpose();
        Dvv.noalias() += mij * ij.transpose();

        Dov.row(i).noalias() += rrij * lev.col(j);
        Dov.row(i).noalias() += iij * exl.col(j);

        virt.noalias() += lij * ij.transpose();
        if (i != j) {
          Dvv.noalias() += lij.transpose() * rij;
          Dvv.noalias() += mij.transpose() * ij;

          Dov.row(j).noalias() += rrij.transpose() * lev.col(i);
          Dov.row(j).noalias() += iij.transpose() * exl.col(i);

          virt.noalias() += lij.transpose() * ij;
        }
      }
    }
    Dov.noalias() -= rev.transpose() * virt;

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
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, 0);
        Doo.noalias() -= rab * lab.transpose();
        Doo.noalias() -= ab * mab.transpose();
        occ.noalias() -= ab * lab.transpose();
        if (a != b) {
          Doo.noalias() -= rab.transpose() * lab;
          Doo.noalias() -= ab.transpose() * mab;
          occ.noalias() -= ab.transpose() * lab;
        }
      }
    }
    Dov.noalias() += occ * rev.transpose();

    // Done with ab-wise amplitudes, must order integrals back.
    this->reorderTensor(*_Jia, _no, _nv);
    if (_Jia != _Bia) {
      this->reorderTensor(*_Bia, _no, _nv);
    }
    this->reorderTensor(_Fia, _no, _nv);
    this->reorderTensor(_gsLagrange, _no, _nv);
  } /* iEigen */

  Timings::timeTaken("CC2 -        Density Matrices");
} /* this->calculateExcitedStateDensities() restricted */

template<>
void CC2Controller<UNRESTRICTED>::calculateExcitedStateDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                 Eigen::VectorXd eigenvalues,
                                                                 std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");
  this->calculateGroundStateDensity();

  unsigned nEigen = _settings.ccexdens ? eigenvalues.size() : 0;
  unsigned mo_a = _no_a + _nv_a;
  unsigned mo_b = _no_b + _nv_b;
  densityMatrices[2] = Eigen::MatrixXd::Zero(mo_a * mo_a + mo_b * mo_b, 1 + nEigen);

  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, nEigen);

  // First column: only ground-state density.
  densityMatrices[2].topRows(mo_a * mo_a).colwise() +=
      Eigen::Map<Eigen::VectorXd>(_gsDensity.alpha.data(), _gsDensity.alpha.size());
  densityMatrices[2].bottomRows(mo_b * mo_b).colwise() +=
      Eigen::Map<Eigen::VectorXd>(_gsDensity.beta.data(), _gsDensity.beta.size());

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[2].col(1 + iEigen).data(), mo_a, mo_a);
    Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[2].col(1 + iEigen).data() + mo_a * mo_a, mo_b, mo_b);

    Eigen::MatrixXd virt_a = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
    Eigen::MatrixXd occ_a = Eigen::MatrixXd::Zero(_no_a, _no_a);

    Eigen::MatrixXd virt_b = Eigen::MatrixXd::Zero(_nv_b, _nv_b);
    Eigen::MatrixXd occ_b = Eigen::MatrixXd::Zero(_no_b, _no_b);

    Eigen::Ref<Eigen::MatrixXd> Doo_a = D_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> Dov_a = D_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> Dvo_a = D_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> Dvv_a = D_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> Doo_b = D_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> Dov_b = D_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> Dvo_b = D_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> Dvv_b = D_b.bottomRightCorner(_nv_b, _nv_b);

    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[2].col(1 + iEigen);
    double eigenvalue = eigenvalues(iEigen);

    Eigen::Map<Eigen::MatrixXd> rev_a(rightEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> lev_a(leftEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> exl_a(excLagrange.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fai_a(_Fai.data(), _nv_a, _no_a);

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

    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
    this->performTransformation(_Yia.alpha, excLagrange.head(_alpha), true, 1);
    this->performTransformation(_Yia.beta, excLagrange.tail(_beta), true, -1);
    this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);

    Eigen::MatrixXd rg_a = -1 * rev_a.transpose() * lev_a;
    Eigen::MatrixXd gr_a = -1 * lev_a * rev_a.transpose();
    Eigen::MatrixXd rg_b = -1 * rev_b.transpose() * lev_b;
    Eigen::MatrixXd gr_b = -1 * lev_b * rev_b.transpose();
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> wia_a(_Wia.alpha.col(Q).data(), _nv_a, _no_a);
      wia_a = jia_a * rg_a + gr_a * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> wia_b(_Wia.beta.col(Q).data(), _nv_b, _no_b);
      wia_b = jia_b * rg_b + gr_b * jia_b;
    }

    _leftSingles = leftEigenvector;
    _exlSingles1 = excLagrange;
    _exlSingles2 = leftEigenvector;

    Doo_a.noalias() -= rev_a.transpose() * lev_a;
    Dvv_a.noalias() += lev_a * rev_a.transpose();
    Dvo_a.noalias() += exl_a;

    Doo_b.noalias() -= rev_b.transpose() * lev_b;
    Dvv_b.noalias() += lev_b * rev_b.transpose();
    Dvo_b.noalias() += exl_b;

    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();

        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, 0, 1);

        Dvv_a.noalias() += lij * rij.transpose();
        Dvv_a.noalias() += mij * ij.transpose();

        Dov_a.row(i).noalias() += rrij * lev_a.col(j);
        Dov_a.row(i).noalias() += iij * exl_a.col(j);

        virt_a.noalias() += lij * ij.transpose();
        if (i != j) {
          Dvv_a.noalias() += lij.transpose() * rij;
          Dvv_a.noalias() += mij.transpose() * ij;

          Dov_a.row(j).noalias() += rrij.transpose() * lev_a.col(i);
          Dov_a.row(j).noalias() += iij.transpose() * exl_a.col(i);

          virt_a.noalias() += lij.transpose() * ij;
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, -1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();

        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, 0, -1);

        Dvv_b.noalias() += lij * rij.transpose();
        Dvv_b.noalias() += mij * ij.transpose();

        Dov_b.row(i).noalias() += rrij * lev_b.col(j);
        Dov_b.row(i).noalias() += iij * exl_b.col(j);

        virt_b.noalias() += lij * ij.transpose();
        if (i != j) {
          Dvv_b.noalias() += lij.transpose() * rij;
          Dvv_b.noalias() += mij.transpose() * ij;

          Dov_b.row(j).noalias() += rrij.transpose() * lev_b.col(i);
          Dov_b.row(j).noalias() += iij.transpose() * exl_b.col(i);

          virt_b.noalias() += lij.transpose() * ij;
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 0);

        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
        Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, 0, 0);

        Dvv_a.noalias() += lij * rij.transpose();
        Dvv_a.noalias() += mij * ij.transpose();

        Dov_a.row(i).noalias() += _oss * rij * lev_b.col(j);
        Dov_a.row(i).noalias() += _oss * ij * exl_b.col(j);

        virt_a.noalias() += lij * ij.transpose();

        Dvv_b.noalias() += lij.transpose() * rij;
        Dvv_b.noalias() += mij.transpose() * ij;

        Dov_b.row(j).noalias() += _oss * rij.transpose() * lev_a.col(i);
        Dov_b.row(j).noalias() += _oss * ij.transpose() * exl_a.col(i);

        virt_b.noalias() += lij.transpose() * ij;
      }
    }
    Dov_a.noalias() -= rev_a.transpose() * virt_a;
    Dov_b.noalias() -= rev_b.transpose() * virt_b;

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
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 1);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, 0, 1);
        Doo_a.noalias() -= rab * lab.transpose();
        Doo_a.noalias() -= ab * mab.transpose();
        occ_a.noalias() -= ab * lab.transpose();
        if (a != b) {
          Doo_a.noalias() -= rab.transpose() * lab;
          Doo_a.noalias() -= ab.transpose() * mab;
          occ_a.noalias() -= ab.transpose() * lab;
        }
      }
    }
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, -1);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, -1);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, -1);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, 0, -1);
        Doo_b.noalias() -= rab * lab.transpose();
        Doo_b.noalias() -= ab * mab.transpose();
        occ_b.noalias() -= ab * lab.transpose();
        if (a != b) {
          Doo_b.noalias() -= rab.transpose() * lab;
          Doo_b.noalias() -= ab.transpose() * mab;
          occ_b.noalias() -= ab.transpose() * lab;
        }
      }
    }
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 0);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, 0);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 0);
        Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, 0, 0);
        Doo_a.noalias() -= rab * lab.transpose();
        Doo_a.noalias() -= ab * mab.transpose();
        occ_a.noalias() -= ab * lab.transpose();
        Doo_b.noalias() -= rab.transpose() * lab;
        Doo_b.noalias() -= ab.transpose() * mab;
        occ_b.noalias() -= ab.transpose() * lab;
      }
    }
    Dov_a.noalias() += occ_a * rev_a.transpose();
    Dov_b.noalias() += occ_b * rev_b.transpose();

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
} /* this->calculateExcitedStateDensities() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity