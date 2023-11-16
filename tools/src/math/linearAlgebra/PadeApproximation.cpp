/**
 * @file   PadeApproximation.cpp
 * @author Johannes Toelle
 *
 * @date 09.10.2020
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
#include "math/linearAlgebra/PadeApproximation.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>
namespace Serenity {

PadeApproximation::PadeApproximation(Eigen::VectorXcd points, Eigen::VectorXcd functionValues)
  : _points(points), _functionValues(functionValues), _numberOfPoints(_functionValues.size()) {
  _padeCoeffs = this->calculatePadeCoeffs();
}

const Eigen::VectorXcd PadeApproximation::calculatePadeCoeffs() {
  Eigen::MatrixXcd g = Eigen::MatrixXcd::Zero(_numberOfPoints, _numberOfPoints);
  g.row(0) = _functionValues;
  for (unsigned int iN = 1; iN < _numberOfPoints; iN++) {
    for (unsigned int z = iN; z < _numberOfPoints; z++) {
      std::complex<double> tmp1 = g(iN - 1, iN - 1) / g(iN - 1, z);
      std::complex<double> tmp2 = g(iN - 1, z) / g(iN - 1, z);
      g(iN, z) = (tmp1 - tmp2) / (_points(z) - _points(iN - 1));
    }
  }
  const Eigen::VectorXcd coeffs = g.diagonal();
  return coeffs;
}

const std::complex<double> PadeApproximation::padeApproximation(std::complex<double> value) {
  Eigen::VectorXcd a = Eigen::VectorXcd::Zero(_numberOfPoints + 1);
  Eigen::VectorXcd b = Eigen::VectorXcd::Zero(_numberOfPoints + 1);
  a(1) = std::complex<double>(_padeCoeffs(0));
  b(0) = std::complex<double>(1.0, 0.0);
  b(1) = std::complex<double>(1.0, 0.0);
  for (unsigned int i = 2; i < _numberOfPoints + 1; i++) {
    a(i) = a(i - 1) + ((value - _points(i - 2)) * _padeCoeffs(i - 1) * a(i - 2));
    b(i) = b(i - 1) + ((value - _points(i - 2)) * _padeCoeffs(i - 1) * b(i - 2));
  }
  std::complex<double> funcValue = a(_numberOfPoints) / b(_numberOfPoints);
  return funcValue;
}

} /* namespace Serenity */