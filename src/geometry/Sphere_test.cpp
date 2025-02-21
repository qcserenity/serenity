/**
 * @file Sphere_test.cpp
 *
 * @date: Feb 28, 2017
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
#include "geometry/Sphere.h"
#include "geometry/Point.h"
#include "geometry/Triangle.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

TEST(SphereTest, simple) {
  Point p1(1, 1, 1);
  Sphere s1(p1, 1.0, SphereType::Ghost, 0);

  Point p2(1, 1, -1);
  Sphere s2(p2, 2.0, SphereType::Ghost, 0);

  auto surf = s1.getSurface();

  double area = 0.0;
  for (auto tri : surf) {
    area += tri.getArea();
  }

  auto surf2 = s2.getSurface();

  double area2 = 0.0;
  for (auto tri : surf2) {
    area2 += tri.getArea();
  }

  unsigned int surfSize = surf.size();
  unsigned int surf2Size = surf2.size();
  unsigned int sixty = 60;

  EXPECT_EQ(sixty, surfSize);
  EXPECT_EQ(sixty, surf2Size);
  EXPECT_NEAR(2.0, s2.distanceTo(s1), 1e-6);
  EXPECT_NEAR(11.396702, area, 1e-6);
  EXPECT_NEAR(45.586808, area2, 1e-6);
}
} // namespace Serenity
