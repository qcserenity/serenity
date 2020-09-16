/**
 * @file   VectorMaths_test.cpp
 *
 * @date   Oct 2, 2014
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
#include "math/VectorMaths.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests VectorMaths.h: addition of vectors
 */
TEST(VectorMathsTest, AdditionOneWorks) {
  std::vector<double> vecOne(3, 2.0);
  std::vector<double> vecTwo(3, 1.0);
  auto result = vecOne + vecTwo;
  EXPECT_EQ(result[0], 3.0);
  EXPECT_EQ(result[1], 3.0);
  EXPECT_EQ(result[2], 3.0);
  EXPECT_NE(result[0], 4.0);
}

/**
 * @test
 * @brief Tests VectorMaths.h: addition of vectors
 */
TEST(VectorMathsTest, AdditionTwoWorks) {
  std::vector<double> vecOne(3, 3.0);
  std::vector<double> vecTwo(3, 1.0);

  auto result = addToVector(vecOne, vecTwo);
  EXPECT_EQ(result[0], 4.0);
  EXPECT_EQ(result[1], 4.0);
  EXPECT_EQ(result[2], 4.0);
  EXPECT_NE(result[0], 5.0);
}

/**
 * @test
 * @brief Tests VectorMaths.h: death test addition of vectors
 */
TEST(VectorMathsTest, AdditionOneFails) {
  std::vector<double> vecOne(3, 3.0);
  std::vector<double> vecTwo(4, 1.0);
  EXPECT_THROW(vecOne + vecTwo, SerenityError);
}

/**
 * @test
 * @brief Tests VectorMaths.h: death tests addition of vectors
 */
TEST(VectorMathsTest, AdditionTwoFails) {
  std::vector<double> vecOne(3, 3.0);
  std::vector<double> vecTwo(4, 1.0);
  EXPECT_THROW(addToVector(vecOne, vecTwo), SerenityError);
}

} /* namespace Serenity */
