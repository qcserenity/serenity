/**
 * @file CartesianToSphericalTransformer_test.cpp
 *
 * @date Jun 13, 2018
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

/* Include Serenity Internal Headers */
#include "basis/CartesianToSphericalTransformer.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class CartesianToSphericalTransformerTest : public ::testing::Test {
 protected:
  CartesianToSphericalTransformerTest() {
  }

  virtual ~CartesianToSphericalTransformerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test CartesianToSphericalTransformerTest
 * @brief Tests the CartesianToSphericalTransformer.h/.cpp for l=2;
 */
TEST_F(CartesianToSphericalTransformerTest, AngularMomemtum_d) {
  unsigned int l = 2;
  auto m = CartesianToSphericalTransformer::getTransformationMatrix(l);

  Eigen::MatrixXd d(6, 5);
  d << 0.000000000000e+00, 0.000000000000e+00, -5.000000000000e-01, 0.000000000000e+00, 8.660254037844e-01,
      1.732050807569e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 1.732050807569e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, -5.000000000000e-01, 0.000000000000e+00, -8.660254037844e-01,
      0.000000000000e+00, 1.732050807569e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, 1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00;

  double test = (d - m).sum();
  EXPECT_NEAR(test, 0.0, 1e-5);
}

/**
 * @test CartesianToSphericalTransformerTest
 * @brief Tests the CartesianToSphericalTransformer.h/.cpp for l=3;
 */
TEST_F(CartesianToSphericalTransformerTest, AngularMomemtum_f) {
  unsigned int l = 3;
  auto m = CartesianToSphericalTransformer::getTransformationMatrix(l);
  Eigen::MatrixXd d(10, 7);
  d << 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, -6.123724356958e-01,
      0.000000000000e+00, 7.905694150421e-01, 2.371708245126e+00, 0.000000000000e+00, -6.123724356958e-01,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, -1.500000000000e+00, 0.000000000000e+00, 1.936491673104e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      -6.123724356958e-01, 0.000000000000e+00, -2.371708245126e+00, 0.000000000000e+00, 3.872983346207e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 2.449489742783e+00,
      0.000000000000e+00, 0.000000000000e+00, -7.905694150421e-01, 0.000000000000e+00, -6.123724356958e-01,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, -1.500000000000e+00, 0.000000000000e+00, -1.936491673104e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 2.449489742783e+00, 0.000000000000e+00,
      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
      0.000000000000e+00, 1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00;
  double test = (d - m).sum();
  EXPECT_NEAR(test, 0.0, 1e-5);
}

} /* namespace Serenity */
