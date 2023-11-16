/**
 * @file CC2DensityGroundstate.cpp
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
void CC2Controller<RESTRICTED>::calculateGroundStateDensity() {
  // SCF occupations.
  for (size_t i = 0; i < _no; ++i) {
    _gsDensity(i, i) = 1.0;
  }

  Eigen::Ref<Eigen::MatrixXd> Doo = _gsDensity.topLeftCorner(_no, _no);
  Eigen::Ref<Eigen::MatrixXd> Dov = _gsDensity.topRightCorner(_no, _nv);
  Eigen::Ref<Eigen::MatrixXd> Dvo = _gsDensity.bottomLeftCorner(_nv, _no);
  Eigen::Ref<Eigen::MatrixXd> Dvv = _gsDensity.bottomRightCorner(_nv, _nv);

  Dvo += Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data(), _nv, _no);

  this->performTransformation(_Zia, _gsLagrange, true);
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      auto tij = this->getGLagrangeAmplitudes(i, j);
      auto ij = this->getAmplitudesA(i, j);
      Eigen::MatrixXd iij = _soss * ij - _sss * ij.transpose();
      Dvv.noalias() += tij * ij.transpose();
      Dov.row(i).noalias() += iij * _gsLagrange.segment(j * _nv, _nv);
      if (i != j) {
        Dvv.noalias() += tij.transpose() * ij;
        Dov.row(j).noalias() += iij.transpose() * _gsLagrange.segment(i * _nv, _nv);
      }
    }
  }

  this->reorderTensor(_Fia, _nv, _no);
  this->reorderTensor(_gsLagrange, _nv, _no);
  this->reorderTensor(_Zia, _nv, _no);
  this->reorderTensor(*_Bia, _nv, _no);
  if (_Jia != _Bia) {
    this->reorderTensor(*_Jia, _nv, _no);
  }

  for (size_t a = 0; a < _nv; ++a) {
    for (size_t b = a; b < _nv; ++b) {
      auto tab = this->getGLagrangeAmplitudesV(a, b);
      auto ab = this->getAmplitudesAV(a, b);
      Doo.noalias() -= ab * tab.transpose();
      if (a != b) {
        Doo.noalias() -= ab.transpose() * tab;
      }
    }
  }

  this->reorderTensor(_Fia, _no, _nv);
  this->reorderTensor(_gsLagrange, _no, _nv);
  this->reorderTensor(*_Bia, _no, _nv);
  if (_Jia != _Bia) {
    this->reorderTensor(*_Jia, _no, _nv);
  }

  // Restricted SCF.
  _gsDensity *= 2.0;

  Eigen::EigenSolver<Eigen::MatrixXd> es(_gsDensity);
  Eigen::VectorXd occupations = es.eigenvalues().real();

  // Sort in descending order and print occupation numbers.
  printf("\n");
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printTableHead("  CC2 unrelaxed natural occ. numbers (1e-2+):");
  }
  else {
    printTableHead("  MP2 unrelaxed natural occ. numbers (1e-2+):");
  }
  unsigned iMax;
  for (unsigned i = 0; i < occupations.size(); ++i) {
    occupations.tail(occupations.size() - i).maxCoeff(&iMax);
    occupations.row(i).swap(occupations.row(iMax + i));
    const char* str = (i % 5 == 0) ? "\n   %8.4f" : "   %8.4f";
    if (occupations(i) > 1e-2) {
      printf(str, occupations(i));
    }
  }
  printf("\n\n  Sum of all occupations %8.4f\n\n", occupations.sum());
}

template<>
void CC2Controller<UNRESTRICTED>::calculateGroundStateDensity() {
  // SCF occupations.
  for (size_t i = 0; i < _no_a; ++i) {
    _gsDensity.alpha(i, i) = 1.0;
  }
  for (size_t i = 0; i < _no_b; ++i) {
    _gsDensity.beta(i, i) = 1.0;
  }

  Eigen::Ref<Eigen::MatrixXd> Doo_a = _gsDensity.alpha.topLeftCorner(_no_a, _no_a);
  Eigen::Ref<Eigen::MatrixXd> Dov_a = _gsDensity.alpha.topRightCorner(_no_a, _nv_a);
  Eigen::Ref<Eigen::MatrixXd> Dvo_a = _gsDensity.alpha.bottomLeftCorner(_nv_a, _no_a);
  Eigen::Ref<Eigen::MatrixXd> Dvv_a = _gsDensity.alpha.bottomRightCorner(_nv_a, _nv_a);

  Eigen::Ref<Eigen::MatrixXd> Doo_b = _gsDensity.beta.topLeftCorner(_no_b, _no_b);
  Eigen::Ref<Eigen::MatrixXd> Dov_b = _gsDensity.beta.topRightCorner(_no_b, _nv_b);
  Eigen::Ref<Eigen::MatrixXd> Dvo_b = _gsDensity.beta.bottomLeftCorner(_nv_b, _no_b);
  Eigen::Ref<Eigen::MatrixXd> Dvv_b = _gsDensity.beta.bottomRightCorner(_nv_b, _nv_b);

  Dvo_a = Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data(), _nv_a, _no_a);
  Dvo_b = Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data() + _alpha, _nv_b, _no_b);

  this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
  this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);

  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      auto tij = this->getGLagrangeAmplitudes(i, j, 1);
      auto ij = this->getAmplitudesA(i, j, 1);
      Eigen::MatrixXd ij_ = _sss * ij - _sss * ij.transpose();
      Dvv_a.noalias() += tij * ij.transpose();
      Dov_a.row(i).noalias() += ij_ * _gsLagrange.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        Dvv_a.noalias() += tij.transpose() * ij;
        Dov_a.row(j).noalias() += ij_.transpose() * _gsLagrange.middleRows(i * _nv_a, _nv_a);
      }
    }
  }

  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      auto tij = this->getGLagrangeAmplitudes(i, j, -1);
      auto ij = this->getAmplitudesA(i, j, -1);
      Eigen::MatrixXd ij_ = _sss * ij - _sss * ij.transpose();
      Dvv_b.noalias() += tij * ij.transpose();
      Dov_b.row(i).noalias() += ij_ * _gsLagrange.middleRows(_alpha + j * _nv_b, _nv_b);
      if (i != j) {
        Dvv_b.noalias() += tij.transpose() * ij;
        Dov_b.row(j).noalias() += ij_.transpose() * _gsLagrange.middleRows(_alpha + i * _nv_b, _nv_b);
      }
    }
  }

  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      auto tij = this->getGLagrangeAmplitudes(i, j, 0);
      auto ij = this->getAmplitudesA(i, j, 0);
      Dvv_a.noalias() += tij * ij.transpose();
      Dov_a.row(i).noalias() += _oss * ij * _gsLagrange.middleRows(_alpha + j * _nv_b, _nv_b);
      Dvv_b.noalias() += tij.transpose() * ij;
      Dov_b.row(j).noalias() += _oss * ij.transpose() * _gsLagrange.middleRows(i * _nv_a, _nv_a);
    }
  }

  this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
  this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
  this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
  this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);
  this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
  this->reorderTensor(_Zia.beta, _nv_b, _no_b);
  this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
  this->reorderTensor(_Bia->beta, _nv_b, _no_b);
  if (_Jia != _Bia) {
    this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
    this->reorderTensor(_Jia->beta, _nv_b, _no_b);
  }

  for (size_t a = 0; a < _nv_a; ++a) {
    for (size_t b = a; b < _nv_a; ++b) {
      auto tab = this->getGLagrangeAmplitudesV(a, b, 1);
      auto ab = this->getAmplitudesAV(a, b);
      Doo_a.noalias() -= ab * tab.transpose();
      if (a != b) {
        Doo_a.noalias() -= ab.transpose() * tab;
      }
    }
  }
  for (size_t a = 0; a < _nv_b; ++a) {
    for (size_t b = a; b < _nv_b; ++b) {
      auto tab = this->getGLagrangeAmplitudesV(a, b, -1);
      auto ab = this->getAmplitudesAV(a, b, -1);
      Doo_b.noalias() -= ab * tab.transpose();
      if (a != b) {
        Doo_b.noalias() -= ab.transpose() * tab;
      }
    }
  }
  for (size_t a = 0; a < _nv_a; ++a) {
    for (size_t b = 0; b < _nv_b; ++b) {
      auto tab = this->getGLagrangeAmplitudesV(a, b, 0);
      auto ab = this->getAmplitudesAV(a, b, 0);
      Doo_a.noalias() -= ab * tab.transpose();
      Doo_b.noalias() -= ab.transpose() * tab;
    }
  }

  this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
  this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
  this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
  this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
  this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
  this->reorderTensor(_Bia->beta, _no_b, _nv_b);
  if (_Jia != _Bia) {
    this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
    this->reorderTensor(_Jia->beta, _no_b, _nv_b);
  }

  Eigen::EigenSolver<Eigen::MatrixXd> es_a(_gsDensity.alpha);
  Eigen::VectorXd occupations = es_a.eigenvalues().real();

  // Sort in descending order and print occupation numbers.
  printf("\n");
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printTableHead("  CC2 unrelaxed natural occ. numbers (1e-2+):");
  }
  else {
    printTableHead("  MP2 unrelaxed natural occ. numbers (1e-2+):");
  }
  unsigned iMax;
  // Sort in descending order and print occupation numbers.
  for (unsigned i = 0; i < occupations.size(); ++i) {
    occupations.tail(occupations.size() - i).maxCoeff(&iMax);
    occupations.row(i).swap(occupations.row(iMax + i));
    const char* str = (i % 5 == 0) ? "\n   %8.4f" : "   %8.4f";
    if (occupations(i) > 1e-2) {
      printf(str, occupations(i));
    }
  }
  printf("\n\n  Sum of all occupations %8.4f\n\n", occupations.sum());

  Eigen::EigenSolver<Eigen::MatrixXd> es_b(_gsDensity.beta);
  occupations = es_b.eigenvalues().real();

  for (unsigned i = 0; i < occupations.size(); ++i) {
    occupations.tail(occupations.size() - i).maxCoeff(&iMax);
    occupations.row(i).swap(occupations.row(iMax + i));
    const char* str = (i % 5 == 0) ? "\n   %8.4f" : "   %8.4f";
    if (occupations(i) > 1e-2) {
      printf(str, occupations(i));
    }
  }
  printf("\n\n  Sum of all occupations %8.4f\n\n", occupations.sum());
}

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity