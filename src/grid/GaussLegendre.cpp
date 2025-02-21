/**
 * @file   GaussLegendre.cpp
 * @author Johannes Toelle
 *
 * @date 09.04.2020
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
#include "grid/GaussLegendre.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>

namespace Serenity {

GaussLegendre::GaussLegendre(unsigned int numberOfPoints) : _numberOfPoints(numberOfPoints) {
  Timings::takeTime("Tech. - Legendre roots/weights");
  // Newton-Raphson convergence threshold is set to numerical limit
  this->legendreWeightsandRoots(_numberOfPoints, 1e-16);
  Timings::timeTaken("Tech. - Legendre roots/weights");
}

const Eigen::VectorXd GaussLegendre::getGridPoints() {
  return _roots;
}
const Eigen::VectorXd GaussLegendre::getWeights() {
  return _weights;
}

std::vector<double> GaussLegendre::legendrePolynomial_derivative(unsigned int order, double x) {
  std::vector<double> legendre_dlegendre;
  legendre_dlegendre.resize(2);
  legendre_dlegendre[0] = x;
  legendre_dlegendre[1] = 0.0;

  double value_minus_1 = 1.0;
  const double f = 1.0 / (pow(x, 2) - 1.0);
  for (unsigned int step = 2; step <= order; step++) {
    const double value = ((2.0 * step - 1.0) * x * legendre_dlegendre[0] - (step - 1.0) * value_minus_1) / step;
    legendre_dlegendre[1] = step * f * (x * value - legendre_dlegendre[0]);

    value_minus_1 = legendre_dlegendre[0];
    legendre_dlegendre[0] = value;
  }
  return legendre_dlegendre;
}

void GaussLegendre::legendreWeightsandRoots(unsigned int order, double threshold) {
  if (order < 2) {
    throw SerenityError("You can not compute the roots for legendre polynomial order = 1!");
  }
  _roots.resize(order);
  _roots.setZero();
  _weights.resize(order);
  _weights.setZero();

  for (unsigned int step = 0; step < order; step++) {
    double root = cos(PI * (step - 0.25) / (order + 0.5));
    auto result = legendrePolynomial_derivative(order, root);

    double newtonRaphsonRatio;
    do {
      newtonRaphsonRatio = result[0] / result[1];
      root -= newtonRaphsonRatio;
      result = legendrePolynomial_derivative(order, root);
    } while (fabs(newtonRaphsonRatio) > threshold);

    _roots(step) = root;
    // Set the second root to be negative
    if (step == 1)
      _roots(step) = -1.0 * _roots(step);
    _weights(step) = 2.0 / ((1 - pow(root, 2)) * pow(result[1], 2));
  }
}

} /* namespace Serenity */