/**
 * @file RegularRankFourTensor_test.cpp
 *
 * @date Feb 19, 2014
 * @author Jan Unsleber
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
#include "math/RegularRankFourTensor.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class RegularRankFourTensorTest
 * @brief Sets everything up for the tests of RegularRankFourTensor.h .
 */
class RegularRankFourTensorTest : public ::testing::Test {
 protected:
  RegularRankFourTensorTest() : tRegularRankFourTensor(4, 1.0) {
    tRegularRankFourTensor(1, 1, 1, 1) = 15.0;
  }

  virtual ~RegularRankFourTensorTest() = default;
  /// the object under test
  RegularRankFourTensor<double> tRegularRankFourTensor;
};

/**
 * @test
 * @brief Tests RegularRankFourTensor.h: retrieve value form the tensor
 */
TEST_F(RegularRankFourTensorTest, GetTensorElement) {
  EXPECT_EQ(tRegularRankFourTensor(1, 1, 1, 1), 15.0);
  EXPECT_EQ(tRegularRankFourTensor(3, 3, 3, 3), 1.0);
  EXPECT_NE(tRegularRankFourTensor(3, 3, 3, 3), 10.0);
}

} // namespace Serenity
