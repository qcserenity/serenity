/**
 * @file FloatMaths_test.cpp
 *
 * @date   Feb 18, 2014
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
#include "math/FloatMaths.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests ::isEqual of cases that should be True.
 *        (including positive, negative and 0 as double).
 */
TEST(FloatMathsTest, True) {
  EXPECT_TRUE(isEqual(0.0, 0.0));
  EXPECT_TRUE(isEqual(5.7564, 5.7564));
  EXPECT_TRUE(isEqual(10000.0, 10000.0));
  EXPECT_TRUE(isEqual(-10.012345, -10.012345));
  EXPECT_TRUE(isEqual(0.000987654321, 0.000987654321));
}

/**
 * @test
 * @brief Tests ::isEqual of cases that should be False
 *       (including positive, negative and 0 as double).
 */
TEST(FloatMathsTest, False) {
  EXPECT_FALSE(isEqual(5.0, 0.0));
  EXPECT_FALSE(isEqual(0.0, 5.0));
  EXPECT_FALSE(isEqual(5.0, -5.0));
  EXPECT_FALSE(isEqual(17.0, 18.0));
  EXPECT_FALSE(isEqual(-7.0, -8.0));
  EXPECT_FALSE(isEqual(0.098765, 0.098763));
  EXPECT_FALSE(isEqual(0.000987654321, 0.000987654322));
}

/**
 * @test
 * @brief Tests ::isEqual with a variable precision
 */
TEST(FloatMathsTest, DifferentEpsilon) {
  EXPECT_FALSE(isEqual(0.9875, 0.9877, 1e-4));
  EXPECT_TRUE(isEqual(0.9875, 0.9877, 1e-3));
}

/**
 * @test
 * @brief Tests ::fipow with powers of 0 to 2 with positive and negative bases
 */
TEST(FloatMathsTest, fipow) {
  EXPECT_EQ(1.0, fipow(-5.234, 0));
  EXPECT_EQ(1.0, fipow(3.214, 0));
  EXPECT_EQ(1.44, fipow(1.2, 2));
  EXPECT_EQ(-7.32, fipow(-7.32, 1));
  // TODO adjust expected accuracy
  //  EXPECT_EQ(0.008, fipow(0.2, 3));
}

} /* namespace Serenity */
