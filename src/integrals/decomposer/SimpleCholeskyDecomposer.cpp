/**
 * @file   SimpleCholeskyDecomposer.cpp
 *
 * @date   Mar 2, 2021
 * @author Lars Hellmann
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
#include "integrals/decomposer/SimpleCholeskyDecomposer.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

SimpleCholeskyDecomposer::SimpleCholeskyDecomposer(Eigen::MatrixXd mat, double thresh)
  : _thresh(thresh), _mat(mat), _vec(mat), _diag(mat.diagonal()) {
  _vec.setZero();
  if (_mat.rows() != mat.cols())
    throw SerenityError("SimpleDecomposer: Matrix not quadratic");
}

bool SimpleCholeskyDecomposer::decompose() {
  auto orig = _mat;
  unsigned int index = 0;
  unsigned int counter = 0;
  _cbas.clear();

  while (true) {
    if (_diag.maxCoeff(&index) <= _thresh)
      break;
    if (counter > _vec.cols())
      break;

    _cbas.push_back(index);
    double factor = 1.0 / std::sqrt(_diag[index]);
    _vec.col(counter) = factor * _mat.col(index);
    _mat -= _vec.col(counter) * _vec.col(counter).transpose();
    _mat.col(index).setZero();
    _mat.row(index).setZero();

    _diag = _mat.diagonal();

    counter++;
  }
  _vec.conservativeResize(_vec.rows(), counter);

  if (_diag.minCoeff() <= -1e-7) {
    throw SerenityError("Simple decomposition failed");
    return false;
  }
  return true;
}

Eigen::MatrixXd SimpleCholeskyDecomposer::getVectors() {
  if (_cbas.size() == 0) {
    decompose();
  }
  return _vec;
}

std::vector<unsigned int> SimpleCholeskyDecomposer::getCholeskyBasis() {
  if (_cbas.size() == 0) {
    decompose();
  }
  return _cbas;
}

} /* namespace Serenity */
