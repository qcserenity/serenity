/**
 * @file Line_test.cpp
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
#include "geometry/Line.h" //To be tested.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <cmath>         //sqrt.

namespace Serenity {

/**
 * @test
 * @brief Tests Line.h: Line--Line intersection.
 */
TEST(LineTest, LineLineIntersection) {
  Line lineA({0, 0, 0}, {1, 1, 1});
  Line lineB({1, 0, 0}, {1, 1, 1});    // Parallel to A, non-identical..
  Line lineC({1, 1, 0}, {-1, -1, -1}); // Anti-parallel to A.
  Line lineD({5, 8, 0}, {0, 0, 1});    // Non-intersecting, non-parallel with A.
  Line lineE({1, 1, 1}, {6, 3, 1});    // Intersecting with to A.
  std::shared_ptr<Eigen::Vector3d> interAB = lineA.getIntersection(lineB);
  EXPECT_FALSE(interAB);
  if (interAB)
    std::cout << *interAB << std::endl;
  std::shared_ptr<Eigen::Vector3d> interAC = lineA.getIntersection(lineC);
  EXPECT_FALSE(interAC);
  if (interAC)
    std::cout << *interAC << std::endl;
  std::shared_ptr<Eigen::Vector3d> interAD = lineA.getIntersection(lineD);
  EXPECT_FALSE(interAD);
  if (interAD)
    std::cout << *interAD << std::endl;
  std::shared_ptr<Eigen::Vector3d> interAE = lineA.getIntersection(lineE);
  Eigen::Vector3d ref = {1, 1, 1};
  EXPECT_NEAR((ref - *interAE).array().abs().sum(), 0.0, 1e-12);
  std::shared_ptr<Eigen::Vector3d> interAA = lineA.getIntersection(lineA);
  ref = {0, 0, 0};
  EXPECT_NEAR((ref - *interAA).array().abs().sum(), 0.0, 1e-12);
}
/**
 * @test
 * @brief Tests Line.h: Direction getter/normalization of the direction vector.
 */
TEST(LineTest, directionGetter) {
  Line lineA({0, 0, 0}, {1, 1, 1});
  Eigen::Vector3d ref = Eigen::Vector3d::Constant(1.0 / sqrt(3.0));
  EXPECT_NEAR((ref - lineA.getDirection()).array().abs().sum(), 0.0, 1e-12);
}
/**
 * @test
 * @brief Tests Line.h: Direction getter/normalization of the direction vector.
 */
TEST(LineTest, isOnLine) {
  Line lineA({0, 0, 0}, {1, 1, 1});
  Eigen::Vector3d rA = {2.2, 2.2, 2.2};
  EXPECT_TRUE(lineA.isOnLine(rA));

  Line lineB({0, 0, 0}, {0, 1, 0});
  Eigen::Vector3d rB = {0, 2, 0};
  EXPECT_TRUE(lineB.isOnLine(rB));
  EXPECT_FALSE(lineB.isOnLine(rA));

  Line lineC({5, 3, 1}, {0, 0, 4});
  Eigen::Vector3d rC = {5, 3, -4};
  EXPECT_TRUE(lineC.isOnLine(rC));
  EXPECT_FALSE(lineC.isOnLine(rB));
  EXPECT_FALSE(lineC.isOnLine(rA));
}

} /* namespace Serenity */
