/**
 * @file JacobiRotation_test.cpp
 *
 * @date Mar 10, 2019
 * @author David Schnieders
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

/* Include Serenity Internal Headers */
#include "math/linearAlgebra/JacobiRotation.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>

namespace Serenity {
/**
 * @class JacobiRotationTest
 * @brief Tests Jacobi rotations on two vectors
 */
class JacobiRotationTest : public ::testing::Test {
 protected:
  JacobiRotationTest() {
  }
  virtual ~JacobiRotationTest() = default;
};

TEST_F(JacobiRotationTest, rotation) {
  Eigen::VectorXd x1(2);
  x1 << 1, 0;

  Eigen::VectorXd x2(2);
  x2 << 0, 1;

  // Rotation of 45 degree
  JacobiRotation::rotate(x1, x2, M_PI / 4.0);

  EXPECT_NEAR(x1[0], 1.0 / sqrt(2), 1.0e-10);
  EXPECT_NEAR(x1[1], 1.0 / sqrt(2), 1.0e-10);
  EXPECT_NEAR(x2[0], -1.0 / sqrt(2), 1.0e-10);
  EXPECT_NEAR(x2[1], 1.0 / sqrt(2), 1.0e-10);

  // Rotate back
  JacobiRotation::rotate(x1, x2, -M_PI / 4.0);

  EXPECT_NEAR(x1[0], 1.0, 1.0e-10);
  EXPECT_NEAR(x1[1], 0.0, 1.0e-10);
  EXPECT_NEAR(x2[0], 0.0, 1.0e-10);
  EXPECT_NEAR(x2[1], 1.0, 1.0e-10);
}

} // namespace Serenity
