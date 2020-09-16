/**
 * @file   ScalarProductMatrix.cpp
 *
 * @date   Jun 11, 2014
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
#include "math/ScalarProductMatrix.h"
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <cassert>
#include <cstring>
#include <memory>

namespace Serenity {
using namespace std;

ScalarProductMatrix::ScalarProductMatrix(const unsigned int rows, const unsigned int cols)
  : _rows(rows), _cols(cols), _rawData(0), _nStored(0) {
  _scalarProductMatrix.reset(new Matrix<double>(_nStored, _nStored));
}

ScalarProductMatrix::~ScalarProductMatrix() {
  eraseAllData();
}

void ScalarProductMatrix::putNewData(Matrix<double> newData) {
  // Copy in new raw data
  _rawData.push_back(newData);
  ++_nStored;
  assert(_rows == newData.rows());
  assert(_cols == newData.cols());
  /*
   * Append a row and column of scalar products
   */
  _scalarProductMatrix->conservativeResize(_nStored, _nStored);
  for (unsigned int i = 0; i < _nStored; ++i) {
    double prod = _rawData[i].cwiseProduct(_rawData[_nStored - 1]).sum();
    (*_scalarProductMatrix)(_nStored - 1, i) = prod;
    (*_scalarProductMatrix)(i, _nStored - 1) = prod;
  }
}

void ScalarProductMatrix::eraseAllData() {
  _nStored = 0;
  _rawData.resize(0);
  _scalarProductMatrix.reset(nullptr);
}

} /* namespace Serenity */
