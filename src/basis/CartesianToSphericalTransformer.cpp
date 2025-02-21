/**
 * @file CartesianToSphericalTransformer.cpp
 *
 * @date May 24, 2018
 * @author Moritz Bensberg
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
#include "basis/CartesianToSphericalTransformer.h"
/* Include Std and External Headers */
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>

namespace Serenity {
/*
 * Map for the transformation matrices.
 */
std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> CartesianToSphericalTransformer::_transformationMatrices = {
    {0, nullptr},  {1, nullptr},  {2, nullptr},  {3, nullptr},  {4, nullptr},  {5, nullptr},  {6, nullptr},
    {7, nullptr},  {8, nullptr},  {9, nullptr},  {10, nullptr}, {11, nullptr}, {12, nullptr}, {13, nullptr},
    {14, nullptr}, {15, nullptr}, {16, nullptr}, {17, nullptr}, {18, nullptr}, {19, nullptr}, {20, nullptr}};

Eigen::MatrixXd& CartesianToSphericalTransformer::getTransformationMatrix(unsigned int l) {
  assert(l <= 20 && "Angular momentum is not supported!");
  if (!_transformationMatrices[l]) {
    // Construct new transformation matrix.
    unsigned int nCar = (l + 1) * (l + 2) / 2;
    unsigned int nSpher = 2 * l + 1;
    Eigen::MatrixXd transformationMatrix(Eigen::MatrixXd::Zero(nCar, nSpher)); // Dimensions nCar x nSpher
    // Loop over orientations.
    for (int m = -(int)l; m <= (int)l; ++m) {
      // The absoulte of m.
      unsigned int abs_m = (m >= 0) ? (unsigned int)m : (unsigned int)(-m);
      // Get normalization.
      const double n = norm(l, abs_m);
      // Loop over t.
      // If looping over t*2 the maximum of the loop reduces to 2*[(l-|m|)/2] = l-|m|
      unsigned int ttimes2Max = (l - abs_m);
      for (unsigned int ttimes2 = 0; ttimes2 <= ttimes2Max; ttimes2 += 2) {
        // Loop over u.
        for (unsigned int u = 0; u <= ttimes2 / 2; ++u) {
          // Loop over vtimes2.
          double v_m = (m >= 0) ? 0 : 0.5;
          // If looping over v*2 the maximum of the loop reduces to 2*([|m|/2-v_m)]+v_m) = |m|
          // which is of course always an integer.
          for (unsigned int vtimes2 = v_m * 2; vtimes2 <= abs_m; vtimes2 += 2) {
            double c = coef(l, abs_m, ttimes2 / 2, u, vtimes2);
            transformationMatrix(mapToCartHarmonics(l, abs_m, ttimes2, 2 * u, vtimes2), m + l) += n * c;
          } // for vtimes2
        }   // for u
      }     // for ttimes2
    }       // for m
    // save transformation matrix.
    _transformationMatrices[l] = std::make_shared<Eigen::MatrixXd>(transformationMatrix);
  } // if find
  return *_transformationMatrices[l];
}

/* Worker functions */

// double CartesianToSphericalTransformer::norm(unsigned int l, unsigned int abs_m) {
//  // Implements eq. 6.4.49 p. 215 Helgaker, Jørgensen, Olsen
//  return sqrt(2*factorial(l+abs_m)*factorial(l-abs_m)/((abs_m == 0)? 2 : 1))/(intPow(2,abs_m)*factorial(l));
//}

double CartesianToSphericalTransformer::coef(unsigned int l, unsigned int abs_m, unsigned int t, unsigned int u,
                                             unsigned int vtimes2) {
  // Implements eq. 6.4.48 p. 215 Helgaker, Jørgensen, Olsen
  //    First element (v=v_m): t+(v-v_m) = t+0
  // Second element (v=v_m+1):           = t+1
  // ...
  double pre = (isEven(t + vtimes2 / 2)) ? 1.0 : -1.0; // Prefactor (-1)^(...) in eq. 6.4.48
  double c = pre * pow(0.25, t) * boost::math::binomial_coefficient<double>(l, t) *
             boost::math::binomial_coefficient<double>(l - t, abs_m + t) *
             boost::math::binomial_coefficient<double>(t, u) * boost::math::binomial_coefficient<double>(abs_m, vtimes2);
  return c;
}

unsigned int CartesianToSphericalTransformer::mapToCartHarmonics(unsigned int l, unsigned int abs_m, unsigned int ttimes2,
                                                                 unsigned int utimes2, unsigned int vtimes2) {
  int lx = ttimes2 + abs_m - utimes2 - vtimes2;
  int ly = utimes2 + vtimes2;
  int lz = l - ttimes2 - abs_m;
  assert(lx >= 0);
  assert(ly >= 0);
  assert(lz >= 0);
  /*
   * Matrix order
   * lx^l, ..., ly^l,...lz^l
   */
  unsigned int index = 0;
  for (int x = l; x >= 0; --x) {
    for (int y = l - x; y >= 0; --y) {
      int z = l - x - y;
      if (x == lx && y == ly && z == lz)
        return index;
      index++;
    }
  }
  assert(false && "get your mapping right ...");
  return 0;
}

} /* namespace Serenity */
