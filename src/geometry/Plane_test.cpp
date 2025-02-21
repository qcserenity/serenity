/**
 * @file Plane_test.cpp
 *
 * @date   Jun 2, 2020
 * @author Moritz Bensberg
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
#include "geometry/Plane.h" //To be tested.
#include "geometry/Line.h"  //Line definition.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <cmath>         //sqrt().

namespace Serenity {

/**
 * @test
 * @brief Tests Plane.h: Plane--Plane intersection.
 */
TEST(PlaneTest, PlanePlaneIntersection) {
  Plane planeXY({0, 0, 0}, {0, 0, 1});
  Plane planeZY({0, 0, 0}, {1, 0, 0});
  Plane diagonal({0, 0, 1}, {-1, 1, 0});
  Plane shiftedDiagonal({1, 0, 1}, {-1, 1, 0});
  std::shared_ptr<Line> yAxis = planeXY.calculateIntersection(planeZY);
  EXPECT_NEAR(yAxis->getPoint().norm(), 0.0, 1e-12);
  Eigen::Vector3d ref = {0, 1, 0};
  EXPECT_NEAR((yAxis->getDirection() - ref).norm(), 0.0, 1e-12);
  std::shared_ptr<Line> xyDiagonal = planeXY.calculateIntersection(diagonal);
  EXPECT_NEAR(xyDiagonal->getPoint().norm(), 0.0, 1e-12);
  ref = {1, 1, 0};
  ref *= -1.0 / sqrt(2);
  EXPECT_NEAR((xyDiagonal->getDirection() - ref).norm(), 0.0, 1e-12);
  std::shared_ptr<Line> withinBothPlanes = planeXY.calculateIntersection(planeXY);
  ref = {0, 1, 0};
  EXPECT_NEAR(withinBothPlanes->getPoint().norm(), 0.0, 1e-12);
  // Note: The direction is in principle arbitrary up to the point that it has
  // to be parallel to the (identical) planes.
  EXPECT_NEAR((withinBothPlanes->getDirection() - ref).norm(), 0.0, 1e-12);
  // Parallel non-intersecting planes.
  std::shared_ptr<Line> nonIntersecting = diagonal.calculateIntersection(shiftedDiagonal);
  EXPECT_FALSE(nonIntersecting);
}
/**
 * @test
 * @brief Tests Plane.h: Plane--Line intersection.
 */
TEST(PlaneTest, PlaneLineIntersection) {
  const Eigen::Vector3d r = {3, 2, 1};
  const Eigen::Vector3d n = {-1, 1, 0};
  const Eigen::Vector3d dParallel = {1, 1, 0};
  Plane diagonal(r, n);
  Line intersectingLine(r, {1.6, 0.86, 4.9});
  Line parallelLine({0.4, 0.9, -0.2}, dParallel);
  Line lineInPlane(r + dParallel, dParallel);
  std::shared_ptr<Eigen::Vector3d> s = diagonal.calculateIntersection(intersectingLine);
  EXPECT_NEAR(0.0, (*s - r).norm(), 1e-9);
  s = diagonal.calculateIntersection(parallelLine);
  EXPECT_FALSE(s);
  s = diagonal.calculateIntersection(lineInPlane);
  EXPECT_NEAR(0.0, (*s - (r + dParallel)).norm(), 1e-9);
}

} /* namespace Serenity */
