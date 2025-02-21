/**
 * @file Triangle_test.cpp
 * @date Feb 28, 2017
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
#include "geometry/Triangle.h"
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

TEST(TriangleTest, simple) {
  Point p1(0, 0, 0);
  Point p2(1, 0, 0);
  Point p3(0, 1, 0);

  Point p4(1, 1, 1);

  Triangle t(p1, p2, p3);

  auto center = t.getCenter();

  t.translate(p4);

  EXPECT_NEAR(0.3333, center.getX(), 1e-4);
  EXPECT_NEAR(0.3333, center.getY(), 1e-4);
  EXPECT_NEAR(0.0, center.getZ(), 1e-6);
  EXPECT_NEAR(1.3333, t.getCenter().getX(), 1e-4);
  EXPECT_NEAR(1.3333, t.getCenter().getY(), 1e-4);
  EXPECT_NEAR(1.0, t.getCenter().getZ(), 1e-6);
  EXPECT_NEAR(0.5, t.getArea(), 1e-6);
}
} // namespace Serenity
