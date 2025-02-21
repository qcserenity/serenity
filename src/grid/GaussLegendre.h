/**
 * @file   GaussLegendre.h
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
#ifndef GAUSSLEGENDRE_H
#define GAUSSLEGENDRE_H

/* Include Serenity Internal Headers */
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>

namespace Serenity {
/**
 * @class GaussLegendre GaussLegendre.h
 * @brief A class, which calculates Gauss-Legendre nodes and weights for quadrature
 *
 * The general procedure is described here:
 *     https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
 * Template for the implementation is from:
 *     https://thoughts-on-coding.com/2019/04/25/numerical-methods-in-c-part-2-gauss-legendre-integration/
 */
class GaussLegendre {
 public:
  /**
   * @brief Constructor
   * @param numberOfPoints The number of weights and nodes
   */
  GaussLegendre(unsigned int numberOfPoints);
  /// @brief Default destructor.
  virtual ~GaussLegendre() = default;
  /// @brief A function to return the gridpoints(nodes)
  const Eigen::VectorXd getGridPoints();
  /// @brief A function to return the weights
  const Eigen::VectorXd getWeights();

 private:
  ///@brief The number of points
  unsigned int _numberOfPoints;
  ///@brief The roots of the Legendre polynomials (corresponds to the nodes)
  Eigen::VectorXd _roots;
  ///@brief The weights of the Gauss-Legendre procedure
  Eigen::VectorXd _weights;
  /**
   * @brief Calculates the legendre polynomial derivative
   * @param order The order of the legendre polynomial derivative
   * @param x The value for which it is evaluated
   * @return The derivative
   */
  std::vector<double> legendrePolynomial_derivative(unsigned int order, double x);
  ///@brief Calculates the roots and weights via Newton-Raphson (threshold = convergence criterion)
  void legendreWeightsandRoots(unsigned int order, double threshold);
};

} // namespace Serenity

#endif /* GAUSSLEGENDRE_H */