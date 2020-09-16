/**
 * @file MatrixFunctions_test.cpp
 *
 * @date Apr 18, 2019
 * @author: Moritz Bensberg
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

#ifndef MATRIXFUNCTIONS_H_
#define MATRIXFUNCTIONS_H_
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/MatrixFunctions.h" //To be tested.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class MatrixFunctionsTest : public ::testing::Test {};

TEST_F(MatrixFunctionsTest, sqrt) {
  Eigen::MatrixXd m(2, 2);
  m << 33, 24, 24, 57;
  Eigen::MatrixXd m_sqrt = mSqrt_Sym(m);
  auto diff = m - m_sqrt * m_sqrt;
  EXPECT_NEAR(diff.array().abs().maxCoeff(), 0.0, 1e-12);
}

TEST_F(MatrixFunctionsTest, psuedoInversSqrt) {
  Eigen::MatrixXd m(2, 2);
  m << 33, 24, 24, 57;
  Eigen::MatrixXd m_sqrt = mSqrt_Sym(m);
  Eigen::MatrixXd pseInv_sqrt = pseudoInversSqrt_Sym(m);
  auto diff = pseInv_sqrt * m_sqrt - Eigen::MatrixXd::Identity(2, 2);
  EXPECT_NEAR(diff.array().abs().maxCoeff(), 0.0, 1e-12);
}

} /* namespace Serenity */

#endif /* MATRIXFUNCTIONS_H_ */
