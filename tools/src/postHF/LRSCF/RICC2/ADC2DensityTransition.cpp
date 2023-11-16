/**
 * @file ADC2DensityTransition.cpp
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

template<>
void ADC2Controller<RESTRICTED>::calculateTransitionDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                              Eigen::VectorXd eigenvalues,
                                                              std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo = _no + _nv;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo * mo, nEigen);
  printf("  Calculating %3i xi density matrices.\n\n", nEigen);

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
    this->reorderTensor(*_Jia, _nv, _no);
    this->reorderTensor(_Xia, _nv, _no);

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
    this->reorderTensor(*_Jia, _no, _nv);
  } /* iEigen */

  // copy right into left density matrices
  densityMatrices[1] = densityMatrices[0];

  Timings::timeTaken("CC2 -        Density Matrices");
} /* this->calculateDensityMatrices() restricted */

template<>
void ADC2Controller<UNRESTRICTED>::calculateTransitionDensities(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                Eigen::VectorXd eigenvalues,
                                                                std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("CC2 -        Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo_a = _no_a + _nv_a;
  unsigned mo_b = _no_b + _nv_b;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo_a * mo_a + mo_b * mo_b, nEigen);
  printf("  Calculating %3i xi density matrices.\n\n", nEigen);

  for (size_t iEigen = 0; iEigen < nEigen; ++iEigen) {
    Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[0].col(iEigen).data(), mo_a, mo_a);
    Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[0].col(iEigen).data() + mo_a * mo_a, mo_b, mo_b);

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
    this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
    this->reorderTensor(_Xia.alpha, _nv_a, _no_a);
    this->reorderTensor(_Jia->beta, _nv_b, _no_b);
    this->reorderTensor(_Xia.beta, _nv_b, _no_b);

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
    this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
    this->reorderTensor(_Jia->beta, _no_b, _nv_b);
  } /* iEigen */

  // copy right into left density matrices
  densityMatrices[1] = densityMatrices[0];

  Timings::timeTaken("CC2 -        Density Matrices");
} /* this->calculateDensityMatrices() unrestricted */

template class ADC2Controller<Options::SCF_MODES::RESTRICTED>;
template class ADC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity