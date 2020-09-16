/**
 *  @file IntegerMaths_test.cpp
 *
 *  @date Feb 12, 2014
 *  @author Jan Unsleber
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
#include "math/IntegerMaths.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests factorial of 0.
 */
TEST(IntegerMathTest, Factorial_Zero) {
  EXPECT_EQ((unsigned int)1, factorial(0));
}
/**
 * @test
 * @brief Tests factorial of positive numbers.
 */
TEST(IntegerMathTest, Factorial_Positive) {
  EXPECT_EQ((unsigned int)1, factorial(1));
  EXPECT_EQ((unsigned int)2, factorial(2));
  EXPECT_EQ((unsigned int)6, factorial(3));
  EXPECT_EQ((unsigned int)40320, factorial(8));
}
/**
 * @test
 * @brief Tests double_factorial.
 */
TEST(IntegerMathTest, DoubleFactorial) {
  EXPECT_EQ((unsigned int)1, double_factorial(1));
  EXPECT_EQ((unsigned int)1, double_factorial(0));
  EXPECT_EQ((unsigned int)2, double_factorial(2));
  EXPECT_EQ((unsigned int)3, double_factorial(3));
  EXPECT_EQ((unsigned int)8, double_factorial(4));
  EXPECT_EQ((unsigned int)15, double_factorial(5));
  EXPECT_EQ((unsigned int)945, double_factorial(9));
}
/**
 * @test
 * @brief Tests isEven.
 */
TEST(IntegerMathTest, IsEven) {
  EXPECT_FALSE(isEven(1));
  EXPECT_TRUE(isEven(2));
  EXPECT_FALSE(isEven(901));
  EXPECT_TRUE(isEven(900));
}

} /* namespace Serenity */
