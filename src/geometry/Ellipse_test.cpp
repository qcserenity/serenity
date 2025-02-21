/**
 * @file Ellipse_test.cpp
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
#include "geometry/Ellipse.h" //To be tested.
#include "geometry/Line.h"    //Line construction for intersection.
#include "geometry/Plane.h"   //Return value of getPlane().
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <cmath>         //sqrt().

namespace Serenity {

/**
 * @test
 * @brief Tests Ellipse.h: getter.
 */
TEST(EllipseTest, Getter) {
  const Eigen::Vector3d center = {1, 1, 1};
  const Eigen::Vector3d r1 = {-1, 1, -1};
  const Eigen::Vector3d r2 = {-0.5, 0, 0.5};
  Ellipse ellipse(center, r1, r2);
  const Plane& plane = ellipse.getPlane();
  Eigen::Vector3d ref = {0.5, 1, 0.5};
  ref *= 2.0 / sqrt(6.0);
  EXPECT_NEAR((ref - plane.getNormalVector()).norm(), 0.0, 1e-12);
  EXPECT_NEAR((center - plane.getPoint()).norm(), 0.0, 1e-12);
  double areaRef = M_PI * r1.norm() * r2.norm();
  EXPECT_NEAR(areaRef, ellipse.getArea(), 1e-12);
  EXPECT_NEAR((center - ellipse.getCenter()).norm(), 0.0, 1e-12);
  EXPECT_NEAR((r1 - ellipse.getR1()).norm(), 0.0, 1e-12);
  EXPECT_NEAR((r2 - ellipse.getR2()).norm(), 0.0, 1e-12);
}
/**
 * @test
 * @brief Tests Ellipse.h: getIntersectionPoints.
 */
TEST(EllipseTest, LineIntersection) {
  const Eigen::Vector3d center = {1, 1, 1};
  const Eigen::Vector3d r1 = {-1, 1, -1};
  const Eigen::Vector3d r2 = {-0.5, 0, 0.5};
  Ellipse ellipse(center, r1, r2);
  Line midLine(center, r1);
  std::vector<Eigen::Vector3d> points = ellipse.getIntersectionPoints(midLine);
  EXPECT_EQ(points.size(), 2);
  Eigen::Vector3d ref2 = center - r1;
  Eigen::Vector3d ref1 = center + r1;
  if (points.size() > 1) {
    EXPECT_NEAR((ref1 - points[0]).norm(), 0.0, 1e-12);
    EXPECT_NEAR((ref2 - points[1]).norm(), 0.0, 1e-12);
  }
  // Construct more arbitrary cuts.
  Eigen::Vector3d cut1 = center + cos(0) * r1 + sin(0) * r2;
  Eigen::Vector3d cut2 = center + cos(M_PI / 2) * r1 + sin(M_PI / 2) * r2;
  Eigen::Vector3d d = cut1 - cut2;
  Line cuttingLine(cut2, d);
  std::vector<Eigen::Vector3d> cuts = ellipse.getIntersectionPoints(cuttingLine);
  EXPECT_EQ(cuts.size(), 2);
  if (cuts.size() > 1) {
    EXPECT_NEAR((cut1 - cuts[0]).norm(), 0.0, 1e-12);
    EXPECT_NEAR((cut2 - cuts[1]).norm(), 0.0, 1e-12);
  }
  cut1 = center + cos(5) * r1 + sin(5) * r2;
  cut2 = center + cos(2) * r1 + sin(2) * r2;
  d = cut1 - cut2;
  Line cuttingLine2(cut2, d);
  cuts = ellipse.getIntersectionPoints(cuttingLine2);
  EXPECT_EQ(cuts.size(), 2);
  if (cuts.size() > 1) {
    EXPECT_NEAR((cut1 - cuts[1]).norm(), 0.0, 1e-12);
    EXPECT_NEAR((cut2 - cuts[0]).norm(), 0.0, 1e-12);
  }
  // Line touches the ellipse.
  const Line touchingLine(ref1, r2);
  cuts = ellipse.getIntersectionPoints(touchingLine);
  EXPECT_EQ(cuts.size(), 0);

  // Line does not cut the ellipse
  const Line passingLine({0, 0, 0}, {1, 0, 0});
  cuts = ellipse.getIntersectionPoints(passingLine);
  EXPECT_EQ(cuts.size(), 0);
}

/**
 * @test
 * @brief Tests Ellipse.h: cutWithPlane.
 *        The ellipse is placed in the xy-plane.
 */
TEST(EllipseTest, CutWithPlaneOrigin) {
  const Eigen::Vector3d center = {0, 0, 0};
  const Eigen::Vector3d r1 = {1, 0, 0};
  const Eigen::Vector3d r2 = {0, 1, 0};
  const Eigen::Vector3d origin = {0, 0, 0};
  const Eigen::Vector3d reference = {1, 1, 0};
  Ellipse ellipse(center, r1, r2);
  Plane middleCut(center, r2);
  CutEllipse fragment1;
  CutEllipse fragment2;
  bool cutResult = ellipse.cutWithPlane(middleCut, reference, fragment1, fragment2);
  EXPECT_TRUE(cutResult);
  EXPECT_NEAR(fragment1.centerOfGravity(1), 0.424413, 1e-5);
  EXPECT_NEAR(fragment1.centerOfGravity(1), -fragment2.centerOfGravity(1), 1e-5);
  EXPECT_NEAR(fragment1.centerOfGravity(0), 0.0, 1e-5);
  EXPECT_NEAR(fragment1.centerOfGravity(2), 0.0, 1e-5);
  Plane middleCut2(center, r1);
  cutResult = ellipse.cutWithPlane(middleCut2, origin, fragment1, fragment2);
  EXPECT_TRUE(cutResult);
  EXPECT_NEAR(fragment1.centerOfGravity(0), -fragment2.centerOfGravity(0), 1e-5);
  EXPECT_NEAR(fragment1.centerOfGravity(1), 0.0, 1e-5);
  EXPECT_NEAR(fragment1.centerOfGravity(2), 0.0, 1e-5);

  double angle = M_PI / 4.0;
  Eigen::Vector3d cut1 = center + cos(angle) * r1 + sin(angle) * r2;
  Eigen::Vector3d cut2 = -cut1;
  Eigen::Vector3d n = cos(angle) * r1 + sin(-angle) * r2;
  ;
  Plane diagonalCut(center, n);
  Line cuttingLine(center, cut1);
  cutResult = ellipse.cutWithPlane(diagonalCut, origin, fragment1, fragment2);
  auto intersectionPoints = ellipse.getIntersectionPoints(cuttingLine);
  EXPECT_TRUE(cutResult);
  double diff = (intersectionPoints[0] - cut2).norm();
  EXPECT_NEAR(0.0, diff, 1e-6);
  diff = (intersectionPoints[1] - cut1).norm();
  EXPECT_NEAR(0.0, diff, 1e-6);
  EXPECT_NEAR(fragment1.centerOfGravity(0), -fragment1.centerOfGravity(1), 1e-6);
}
/**
 * @test
 * @brief Tests Ellipse.h: cutWithPlane.
 *        More arbitrary ellipse.
 */
TEST(EllipseTest, CutWithPlane) {
  const Eigen::Vector3d center = {1, 1, 1};
  const Eigen::Vector3d r1 = {-1, 1, -1};
  const Eigen::Vector3d r2 = {-0.5, 0, 0.5};
  const Eigen::Vector3d origin = {0, 0, 0};
  const Eigen::Vector3d xAxis = {1, 0, 0};
  Ellipse ellipse(center, r1, r2);
  Plane middleCut(center, r2);
  CutEllipse fragment1;
  CutEllipse fragment2;
  // Cut the middle along r1.
  bool cutResult = ellipse.cutWithPlane(middleCut, origin, fragment1, fragment2);
  EXPECT_TRUE(cutResult);
  EXPECT_NEAR(fragment1.area, fragment2.area, 1e-12);
  double centerOfGrvityTest = (0.5 * (fragment1.centerOfGravity + fragment2.centerOfGravity) - center).array().abs().sum();
  EXPECT_NEAR(centerOfGrvityTest, 0.0, 1e-12);
  // Cut the middle along r2
  Plane middleCut2(center, r1);
  cutResult = ellipse.cutWithPlane(middleCut2, origin, fragment1, fragment2);
  EXPECT_TRUE(cutResult);
  EXPECT_NEAR(fragment1.area, fragment2.area, 1e-12);
  centerOfGrvityTest = (0.5 * (fragment1.centerOfGravity + fragment2.centerOfGravity) - center).array().abs().sum();
  EXPECT_NEAR(centerOfGrvityTest, 0.0, 1e-12);
  // Passing, on the reference side
  Plane passing(center + 1 * r1, r1);
  cutResult = ellipse.cutWithPlane(passing, origin, fragment1, fragment2);
  EXPECT_FALSE(cutResult);
  EXPECT_NEAR(ellipse.getArea(), fragment1.area, 1e-12);
  centerOfGrvityTest = (fragment1.centerOfGravity - center).array().abs().sum() + fragment2.centerOfGravity.norm();
  EXPECT_NEAR(centerOfGrvityTest, 0.0, 1e-12);
  // Passing, on the other side
  Plane passing2(center - r1, r1);
  cutResult = ellipse.cutWithPlane(passing2, center - 2.0 * r1, fragment1, fragment2);
  EXPECT_FALSE(cutResult);
  EXPECT_NEAR(ellipse.getArea(), fragment2.area, 1e-12);
  centerOfGrvityTest = (fragment2.centerOfGravity - center).array().abs().sum() + fragment1.centerOfGravity.norm();
  EXPECT_NEAR(centerOfGrvityTest, 0.0, 1e-12);
  // Random cut.
  Eigen::Vector3d cut1 = center + cos(5) * r1 + sin(5) * r2;
  Eigen::Vector3d cut2 = center + cos(2) * r1 + sin(2) * r2;
  Eigen::Vector3d d = cut1 - cut2;
  Plane randomCut(cut2, d.cross(xAxis));
  cutResult = ellipse.cutWithPlane(randomCut, origin, fragment1, fragment2);
  EXPECT_TRUE(cutResult);
  // I did not calculate this by hand, but it looks reasonable to me.
  // Both center of gravities are within the ellipse and close to the original center.
  EXPECT_NEAR(fragment2.centerOfGravity(0), 0.57032116852707615, 1e-4);
  EXPECT_NEAR(fragment2.centerOfGravity(1), 1.3618979448087845, 1e-4);
  EXPECT_NEAR(fragment2.centerOfGravity(2), 0.7058829418553545, 1e-4);
  EXPECT_NEAR(fragment1.centerOfGravity(0), 1.5146601153622066, 1e-4);
  EXPECT_NEAR(fragment1.centerOfGravity(1), 0.56652638114342102, 1e-4);
  EXPECT_NEAR(fragment1.centerOfGravity(2), 1.3522871223509521, 1e-4);
  const Plane& ellipsePlane = ellipse.getPlane();
  EXPECT_TRUE(ellipsePlane.pointIsInPlane(fragment1.centerOfGravity));
  EXPECT_TRUE(ellipsePlane.pointIsInPlane(fragment2.centerOfGravity));
  centerOfGrvityTest =
      (1.0 / ellipse.getArea() * (fragment2.area * fragment2.centerOfGravity + fragment1.area * fragment1.centerOfGravity) - center)
          .array()
          .abs()
          .sum();
  EXPECT_NEAR(centerOfGrvityTest, 0.0, 1e-12);
  // Plane and ellipse plane are identical.
  cutResult = ellipse.cutWithPlane(ellipsePlane, origin, fragment1, fragment2);
  EXPECT_FALSE(cutResult);
  EXPECT_NEAR(fragment1.centerOfGravity.norm(), 0.0, 1e-12);
  EXPECT_NEAR(fragment2.centerOfGravity.norm(), 0.0, 1e-12);
  EXPECT_NEAR(fragment1.area, 0.0, 1e-12);
  EXPECT_NEAR(fragment2.area, 0.0, 1e-12);
  // Plane and ellipse are parallel
  Plane shiftedEllipsePlane(2.0 * center, ellipsePlane.getNormalVector());
  cutResult = ellipse.cutWithPlane(shiftedEllipsePlane, origin, fragment1, fragment2);
  EXPECT_FALSE(cutResult);
  // The ellipse is on the side of the origin.
  EXPECT_NEAR((fragment1.centerOfGravity - center).norm(), 0.0, 1e-12);
  EXPECT_NEAR(fragment2.centerOfGravity.norm(), 0.0, 1e-12);
  EXPECT_NEAR(fragment1.area - ellipse.getArea(), 0.0, 1e-12);
  EXPECT_NEAR(fragment2.area, 0.0, 1e-12);
}

} /* namespace Serenity */
