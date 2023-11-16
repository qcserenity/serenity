/**
 * @file CC2Normalization.cpp
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

template<>
Eigen::MatrixXd CC2Controller<RESTRICTED>::normalizeEigenvectors(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                 Eigen::VectorXd eigenvalues) {
  Timings::takeTime("CC2 -           Normalization");
  Eigen::MatrixXd norms;

  if (eigenvectors.size() == 1) {
    printf("  Normalizing eigenvectors such that <~R|R> = 1.\n\n");
    norms.resize(eigenvalues.size(), 1);
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      double crr = eigenvector.dot(eigenvector);

      this->performTransformation(_Xia, eigenvector, false);
      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = _settings.triplet ? 0 : i; j < _no; ++j) {
          Eigen::MatrixXd amp1 = _Xia.middleRows(i * _nv, _nv) * _Bia->middleRows(j * _nv, _nv).transpose();
          Eigen::MatrixXd amp2 = _Bia->middleRows(i * _nv, _nv) * _Xia.middleRows(j * _nv, _nv).transpose();
          Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j) + eigenvalue) - _eVirt;

          if (!_settings.triplet) {
            Eigen::MatrixXd rij = (amp1 + amp2).cwiseQuotient(denom);
            Eigen::MatrixXd rrij = _soss * rij - _sss * rij.transpose();
            crr += (i == j ? 0.5 : 1.0) * rrij.cwiseProduct(rij).sum();
          }
          else {
            Eigen::MatrixXd rij = 2 * amp1.cwiseQuotient(denom);
            Eigen::MatrixXd rrij =
                ((_sss + _oss) * amp1 + (_sss - _oss) * amp2 - _sss * (amp1 + amp2).transpose()).cwiseQuotient(denom);
            crr += 0.5 * (rrij).cwiseProduct(rij).sum();
          }
        }
      }

      norms(iEigen, 0) = crr;
      // eigenvector *= 1 / std::sqrt(crr);
    } /* iEigen */
  }
  else {
    printf("  Normalizing eigenvectors such that <L|R> = <R|R> = 1.\n\n");
    norms.resize(eigenvalues.size(), 2);
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
      Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      double clr = leftEigenvector.dot(rightEigenvector);
      double crr = rightEigenvector.dot(rightEigenvector);
      this->performTransformation(_Xia, rightEigenvector, false);
      this->performTransformation(_Zia, leftEigenvector, true);
      _leftSingles = leftEigenvector;

      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }

      norms(iEigen, 0) = crr;
      norms(iEigen, 1) = clr;
    } /* iEigen */
  }
  Timings::timeTaken("CC2 -           Normalization");

  return norms;
} /* this->normalizeEigenvectors() restricted */

template<>
Eigen::MatrixXd CC2Controller<UNRESTRICTED>::normalizeEigenvectors(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                   Eigen::VectorXd eigenvalues) {
  Timings::takeTime("CC2 -           Normalization");
  Eigen::MatrixXd norms;

  if (eigenvectors.size() == 1) {
    printf("  Normalizing eigenvectors such that <~R|R> = 1.\n\n");
    norms.resize(eigenvalues.size(), 1);
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      // get eigenpair
      Eigen::Ref<Eigen::VectorXd> eigenvector = eigenvectors[0].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      // calculate scaling factor
      double crr = eigenvector.dot(eigenvector);
      this->performTransformation(_Xia.alpha, eigenvector.head(_alpha), false, 1);
      this->performTransformation(_Xia.beta, eigenvector.tail(_beta), false, -1);

      // alpha
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
          crr += (i == j ? 0.5 : 1.0) * (_sss * rij - _sss * rij.transpose()).cwiseProduct(rij).sum();
        }
      }

      // beta
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
          crr += (i == j ? 0.5 : 1.0) * (_sss * rij - _sss * rij.transpose()).cwiseProduct(rij).sum();
        }
      }

      // mixed (double counting for alpha -> beta and beta -> alpha)
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
          crr += (0.5 + 0.5) * (_oss * rij).cwiseProduct(rij).sum();
        }
      }

      norms(iEigen, 0) = crr;
    } /* iEigen */
  }
  else {
    printf("  Normalizing eigenvectors such that <L|R> = <R|R> = 1.\n\n");
    norms.resize(eigenvalues.size(), 2);
    for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
      Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
      Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
      double eigenvalue = eigenvalues(iEigen);

      // calculate scaling factor
      double clr = leftEigenvector.dot(rightEigenvector);
      double crr = rightEigenvector.dot(rightEigenvector);
      this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
      this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
      this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
      this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);
      _leftSingles = leftEigenvector;

      // alpha
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }

      // beta
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
          clr += (i == j ? 0.5 : 1.0) * lij.cwiseProduct(rij).sum();
          crr += (i == j ? 0.5 : 1.0) * rij.cwiseProduct(rij).sum();
        }
      }

      // mixed (double counting for alpha -> beta and beta -> alpha)
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
          Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);
          clr += (0.5 + 0.5) * lij.cwiseProduct(rij).sum();
          crr += (0.5 + 0.5) * rij.cwiseProduct(rij).sum();
        }
      }

      norms(iEigen, 0) = crr;
      norms(iEigen, 1) = clr;
    } /* iEigen */
  }
  Timings::timeTaken("CC2 -           Normalization");

  return norms;
} /* this->normalizeEigenvectors() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity