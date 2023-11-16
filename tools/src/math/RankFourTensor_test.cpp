/**
 * @file RankFourTensor_test.cpp
 *
 * @date  Feb 19, 2014
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
#include "math/RankFourTensor.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class RankFourTensorTest
 * @brief Sets everything up for the tests of RankFourTensor.h .
 */
class RankFourTensorTest : public ::testing::Test {
 protected:
  RankFourTensorTest() : tRankFourTensor(2, 2, 2, 2) {
    tRankFourTensor(1, 1, 1, 1) = 15.0;
  }

  virtual ~RankFourTensorTest() = default;
  /// The RankFourTensor under test
  RankFourTensor<double> tRankFourTensor;
};

/**
 * @test
 * @brief Tests RankFourTensor.h: retrieve value form the tensor
 */
TEST_F(RankFourTensorTest, GetTensorElement) {
  EXPECT_EQ(tRankFourTensor(1, 1, 1, 1), 15.0);
}

/**
 * @test
 * @brief Tests RankFourTensor.h: print entire tensor.
 */
TEST_F(RankFourTensorTest, PrintTensor) {
  tRankFourTensor.print();
}

} // namespace Serenity
