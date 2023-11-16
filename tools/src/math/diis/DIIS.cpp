/**
 * @file   DIIS.cpp
 *
 * @date   Nov 18, 2013
 * @author Thomas Dresselhaus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "math/diis/DIIS.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"
#include "math/Matrix.h"
#include "math/linearAlgebra/MatrixFunctions.h"
/* Include Std and External Headers */
#include <algorithm>
#include <cmath>
#include <memory>

namespace Serenity {

DIIS::DIIS(const unsigned maxStore, const bool diskMode, double scaling)
  : _maxStore(maxStore),
    _errorVectors(maxStore + 1),
    _targetVectors(maxStore + 1),
    _nStored(0),
    _diskMode(diskMode),
    _scaling(scaling),
    _cycle(0) {
  if (diskMode) {
    _errorDiskVectors.resize(maxStore + 1);
    _targetDiskVectors.resize(maxStore + 1);
  }
}

unsigned int DIIS::getNVectorsStored() {
  return _nStored;
}

void DIIS::store(const Eigen::Ref<const Eigen::VectorXd>& targetVector, const Eigen::Ref<const Eigen::VectorXd>& newErrorVector) {
  _targetVectors[_nStored] = std::make_unique<Eigen::VectorXd>(targetVector);
  _errorVectors[_nStored] = std::make_unique<Eigen::VectorXd>(newErrorVector);
  ++_nStored;

  if (_nStored > _maxStore) {
    this->shiftVectors();
  }
  ++_cycle;
}
void DIIS::storeMatrix(const Eigen::MatrixXd& parameters, const Eigen::MatrixXd& gradients) {
  this->store(Eigen::Map<const Eigen::VectorXd>(parameters.data(), parameters.cols() * parameters.rows()),
              Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.cols() * gradients.rows()));
}

void DIIS::store(const std::vector<double>& parameters, const std::vector<double>& gradients) {
  this->store(Eigen::Map<const Eigen::VectorXd>(parameters.data(), parameters.size()),
              Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.size()));
}
void DIIS::store(VectorOnDiskStorageController& targetVector, VectorOnDiskStorageController& newErrorVector) {
  auto nameTag = "DIISStored_" + std::to_string(_cycle);
  _targetDiskVectors[_nStored] =
      std::make_unique<VectorOnDiskStorageController>(targetVector, nameTag + targetVector.getHDF5FileName());
  _errorDiskVectors[_nStored] =
      std::make_unique<VectorOnDiskStorageController>(newErrorVector, nameTag + newErrorVector.getHDF5FileName());
  ++_nStored;

  if (_nStored > _maxStore) {
    this->shiftVectors();
  }
  ++_cycle;
}

void DIIS::optimize(Eigen::Ref<Eigen::VectorXd> targetVector, const Eigen::Ref<const Eigen::VectorXd>& newErrorVector) {
  this->store(targetVector, newErrorVector);

  this->initNewB();

  // setup B matrix
  for (unsigned i = 0; i < _nStored; ++i) {
    for (unsigned j = 0; j <= i; ++j) {
      const double prod = (*_errorVectors[i]).dot((*_errorVectors[j]));
      _B(i + 1, j + 1) = prod;
      _B(j + 1, i + 1) = prod;
    }
  }

  _nExcluded = 0;
  while (this->checkConditioning() && _nExcluded < _nStored) {
    this->excludeEquation();
  }

  // solve linear system
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(_B.rows());
  rhs(0) = -1.0;

  Eigen::VectorXd coefficients = _B.colPivHouseholderQr().solve(rhs);

  // extrapolate target vector
  targetVector.setZero();
  for (unsigned i = _nExcluded; i < _nStored; ++i) {
    targetVector += (*_targetVectors[i]) * coefficients[i + 1 - _nExcluded];
  }
}

void DIIS::optimize(VectorOnDiskStorageController& targetVector, VectorOnDiskStorageController& newErrorVector) {
  this->store(targetVector, newErrorVector);

  this->initNewB();

  // setup B matrix
  for (unsigned i = 0; i < _nStored; ++i) {
    auto& tmpI = (*_errorDiskVectors[i]);
    for (unsigned j = 0; j <= i; ++j) {
      auto& tmpJ = (*_errorDiskVectors[j]);
      const double prod = tmpI * tmpJ;
      _B(i + 1, j + 1) = prod;
      _B(j + 1, i + 1) = prod;
    }
  }

  _nExcluded = 0;
  while (this->checkConditioning() && _nExcluded < _nStored) {
    this->excludeEquation();
  }

  // solve linear system
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(_B.rows());
  rhs(0) = -1.0;
  Eigen::VectorXd coefficients = _B.colPivHouseholderQr().solve(rhs);

  // extrapolate target vector
  const auto labelList = targetVector.getLabelList();
  for (const auto& label : labelList) {
    auto optSegment = std::make_shared<Eigen::VectorXd>(*targetVector.getVectorSegment(label));
    optSegment->setZero();
    for (unsigned i = _nExcluded; i < _nStored; ++i) {
      *optSegment += (*_targetDiskVectors[i]->getVectorSegment(label)) * coefficients[i + 1 - _nExcluded];
    }
    targetVector.storeVectorSegment(optSegment, label);
  }
}

void DIIS::reinit() {
  OutputControl::vOut << " ***** Re-initialize DIIS Procedure *****" << std::endl;

  this->cleanUp();
  _B.resize(0, 0);
  _targetVectors.resize(_maxStore + 1);
  _targetDiskVectors.resize(_maxStore + 1);
  _errorVectors.resize(_maxStore + 1);
  _errorDiskVectors.resize(_maxStore + 1);
}

void DIIS::shiftVectors() {
  --_nStored;

  for (unsigned i = 0; i < _nStored; i++) {
    if (_diskMode) {
      _errorDiskVectors[i] = std::move(_errorDiskVectors[i + 1]);
      _targetDiskVectors[i] = std::move(_targetDiskVectors[i + 1]);
    }
    else {
      _errorVectors[i] = std::move(_errorVectors[i + 1]);
      _targetVectors[i] = std::move(_targetVectors[i + 1]);
    }
  }
}

bool DIIS::checkConditioning() {
  // this vastly improves numerical stability
  unsigned block = _B.rows() - 1;
  double norm = _B.bottomRightCorner(block, block).cwiseAbs().maxCoeff();
  if (norm > 1e-12) {
    _B.bottomRightCorner(block, block) *= 1.0 / norm;
  }

  // this prevents large coefficients
  _B.diagonal() *= _scaling;

  // check whether the 'normalized' matrix is invertible
  Eigen::FullPivLU<Eigen::MatrixXd> lu(_B);
  return !lu.isInvertible();
}

void DIIS::excludeEquation() {
  ++_nExcluded;
  unsigned newDim = _nStored - _nExcluded;
  Eigen::MatrixXd blockToKeep = _B.bottomRightCorner(newDim, newDim);
  _B.block(1, 1, newDim, newDim) = blockToKeep;
  _B.conservativeResize(newDim + 1, newDim + 1);
}

void DIIS::initNewB() {
  _B = Eigen::MatrixXd::Zero(_nStored + 1, _nStored + 1);
  _B.col(0).setConstant(-1.0);
  _B.row(0).setConstant(-1.0);
  _B(0, 0) = 0;
}

void DIIS::cleanUp() {
  _B.resize(0, 0);
  if (_nStored > 0) {
    _targetVectors.resize(0);
    _targetDiskVectors.resize(0);
    _errorVectors.resize(0);
    _errorDiskVectors.resize(0);
    _nStored = 0;
  }
}

} /* namespace Serenity */
