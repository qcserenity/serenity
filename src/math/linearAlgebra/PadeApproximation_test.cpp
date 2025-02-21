/**
 * @file   PadeApproximation_test.cpp
 *
 * @date   12.10.2020
 * @author J. Toelle
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

/* Include Serenity Internal Headers */
#include "math/linearAlgebra/PadeApproximation.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

namespace Serenity {
TEST(PadeApproximation, fit) {
  Eigen::VectorXcd points(7);
  points << 1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0;
  Eigen::VectorXd functionValues(points.size());
  for (unsigned int i = 0; i < points.size(); i++) {
    double point = points(i).real();
    functionValues(i) = exp(point);
  }
  auto padeApprox = PadeApproximation(points, functionValues);
  const auto testValue = padeApprox.padeApproximation(std::complex<double>(1.4, 0.0));
  EXPECT_NEAR(testValue.real(), exp(1.4), 1e-8);
}
} /* namespace Serenity */