/**
 * @file   Point_test.cpp
 *
 * @date   Oct 1, 2014
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
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests Point.h: creation and getters of stored data
 */
TEST(PointTest, CreatePointAndGetCoordinates) {
  Point point(12.3, -4.0, 0.0);
  EXPECT_EQ(point.getX(), 12.3);
  EXPECT_EQ(point.getY(), -4.0);
  EXPECT_EQ(point.getZ(), 0.0);
  EXPECT_NE(point.getX(), 0.0);
}

/**
 * @test
 * @brief Tests Point.h: distanceToOrigin()
 */
TEST(PointTest, DistanceToOrigin) {
  Point point(0.0, 3.0, 4.0);
  EXPECT_EQ(point.distanceToOrigin(), 5.0);
}

/**
 * @test
 * @brief Tests Point.h: addition and subtraction of points
 */
TEST(PointTest, AddAndSubtract) {
  Point pointOne(12.3, -4.0, 0.0);
  Point pointTwo(0.0, 3.0, 4.0);
  auto result = pointOne + pointTwo;
  EXPECT_EQ(result.getX(), 12.3);
  EXPECT_EQ(result.getY(), -1.0);
  EXPECT_EQ(result.getZ(), 4.0);
  EXPECT_NE(result.getZ(), 0.0);
  result = result - pointTwo;
  EXPECT_EQ(result.getX(), pointOne.getX());
  EXPECT_EQ(result.getY(), pointOne.getY());
  EXPECT_EQ(result.getZ(), pointOne.getZ());
  EXPECT_NE(result.getZ(), 4.0);
}

/**
 * @test
 * @brief Tests Point.h: distance between points
 */
TEST(PointTest, Distance) {
  Point pointOne(45.0, -4.0, 4.0);
  Point pointTwo(41.0, -1.0, 4.0);
  EXPECT_EQ(distance(pointOne, pointTwo), 5.0);
  EXPECT_NE(distance(pointOne, pointTwo), 0.0);
}

} // namespace Serenity
