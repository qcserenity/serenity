/**
 * @file   SimpleCholeskyDecomposer_test.cpp
 *
 * @date   Apr 04, 2020
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

/* Include Serenity Internal Headers */
#include "integrals/decomposer/SimpleCholeskyDecomposer.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class SimpleCholeskyDecomposerTest : public ::testing::Test {
 protected:
  virtual ~SimpleCholeskyDecomposerTest() = default;

  virtual void TearDown() {
  }
};

TEST_F(SimpleCholeskyDecomposerTest, Reconstruction) {
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(100, 100);
  Eigen::VectorXd diag = Eigen::VectorXd::Ones(100);
  Eigen::MatrixXd posDef = mat * mat.transpose();
  posDef += diag.asDiagonal();

  SimpleCholeskyDecomposer sd(posDef, 0.0);
  auto vec = sd.getVectors();
  sd.getCholeskyBasis();

  Eigen::MatrixXd diff = posDef - vec * vec.transpose();

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-12);
}

} // namespace Serenity
