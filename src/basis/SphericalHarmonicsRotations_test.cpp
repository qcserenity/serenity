/**
 * @file SphericalHarmonicsRotations_test.cpp
 *
 * @author Moritz Bensberg
 * @date Feb 12, 2020
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
#include "basis/SphericalHarmonicsRotations.h" //To be tested.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework
#include <math.h>        //sqrt

namespace Serenity {

class SphericalHarmonicsRotationsTest : public ::testing::Test {
 protected:
  SphericalHarmonicsRotationsTest() {
  }

  virtual ~SphericalHarmonicsRotationsTest() = default;
};

TEST_F(SphericalHarmonicsRotationsTest, PrimitiveProperties) {
  unsigned int l = 3;
  Eigen::MatrixXd m = SphericalHarmonicsRotations::getTransformationMatrix(l, 0, 0, 0);
  // If all Euler angles are zero, the transformation has to be the identity!
  double test = (m - Eigen::MatrixXd::Identity(2 * l + 1, 2 * l + 1)).array().abs().sum();
  EXPECT_NEAR(test, 0.0, 1e-12);
  // Check the J matrix for l=3. This is taken from the original reference
  // Journal of Physics A: Mathematical and Theoretical 40, 1597 (2007).
  Eigen::MatrixXd j3(7, 7);
  j3 << 0, 0, 0, sqrt(10.0) / 4, 0, -sqrt(6.0) / 4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, sqrt(6.0) / 4, 0, sqrt(10.0) / 4,
      0, sqrt(10.0) / 4, 0, sqrt(6.0) / 4, 0, 0, 0, 0, 0, 0, 0, 0, -1.0 / 4, 0, -sqrt(15) / 4, -sqrt(6.0) / 4, 0,
      sqrt(10.0) / 4, 0, 0, 0, 0, 0, 0, 0, 0, -sqrt(15.0) / 4, 0, 1.0 / 4;
  Eigen::MatrixXd j3OnTheFly = SphericalHarmonicsRotations::getJMatrix(3);
  test = (j3OnTheFly - j3).array().abs().sum();
  EXPECT_NEAR(test, 0.0, 1e-12);
}
TEST_F(SphericalHarmonicsRotationsTest, PiHalfPRotation) {
  double piHalf = M_PI / 2.0;
  Eigen::MatrixXd ref(3, 3);
  Eigen::MatrixXd m = SphericalHarmonicsRotations::getTransformationMatrix(1, -piHalf, 0, 0);
  //     px  py pz
  ref << 0, 0, -1, // px -->  pz
      0, 1, 0,     // py -->  py
      1, 0, 0;     // pz --> -px
  // Rotation around the z-axis by pi/2 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
  m = SphericalHarmonicsRotations::getTransformationMatrix(1, 0, piHalf, 0);
  ref << 1, 0, 0, 0, 0, -1, 0, 1, 0;
  // Rotation around the y-axis by -pi/2 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
  m = SphericalHarmonicsRotations::getTransformationMatrix(1, 0, 0, piHalf);
  ref << 0, 0, 1, 0, 1, 0, -1, 0, 0;
  // Rotation around the z-axis by -pi/2 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
}
TEST_F(SphericalHarmonicsRotationsTest, PiFourthPRotation) {
  double anlge = M_PI / 4.0;
  Eigen::MatrixXd ref(3, 3);
  Eigen::MatrixXd m = SphericalHarmonicsRotations::getTransformationMatrix(1, -anlge, 0, 0);
  double oneOverSqrt2 = 1.0 / sqrt(2.0);
  ref << oneOverSqrt2, 0, -oneOverSqrt2, 0, 1, 0, oneOverSqrt2, 0, oneOverSqrt2;
  // Rotation around the y-axis by -pi/4 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
  m = SphericalHarmonicsRotations::getTransformationMatrix(1, 0, anlge, 0);
  ref << 1, 0, 0, 0, oneOverSqrt2, -oneOverSqrt2, 0, oneOverSqrt2, oneOverSqrt2;
  // Rotation around the y-axis by pi/4 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
  m = SphericalHarmonicsRotations::getTransformationMatrix(1, 0, 0, anlge);
  ref << oneOverSqrt2, 0, oneOverSqrt2, 0, 1, 0, -oneOverSqrt2, 0, oneOverSqrt2;
  // Rotation around the z-axis by -pi/4 in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
}
TEST_F(SphericalHarmonicsRotationsTest, PiDRotation) {
  double anlge = M_PI;
  Eigen::MatrixXd ref(5, 5);
  Eigen::MatrixXd m = SphericalHarmonicsRotations::getTransformationMatrix(2, anlge, 0, 0);
  ref << 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1;
  // Rotation around the z-axis by pi in a right hand coordinate system.
  EXPECT_NEAR((ref - m).array().abs().sum(), 0.0, 1e-12);
}
TEST_F(SphericalHarmonicsRotationsTest, Zero19Rotation) {
  double anlge = 0;
  Eigen::MatrixXd m = SphericalHarmonicsRotations::getTransformationMatrix(19, anlge, 0, 0);
  EXPECT_NEAR((Eigen::MatrixXd::Identity(2 * 19 + 1, 2 * 19 + 1) - m).array().abs().sum(), 0.0, 1e-12);
}
} /* namespace Serenity */
