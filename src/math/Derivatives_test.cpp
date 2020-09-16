/**
 * @file Derivatives_test.cpp
 *
 * @date August 28, 2015
 * @author Thomas Dresselhaus
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
#include "math/Derivatives.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <vector>

namespace Serenity {

class DerivativesTest : public ::testing::Test {
 protected:
  DerivativesTest()
    : tGradDouble{-1.0, 2.4, 0.1},
      tHessDouble{1.2, 3.4, -4.0, 0.2, 0.5, 0.7},
      tGradIntVec{{1, 2}, {-1, -1}, {5, 7}},
      tHessIntVec{{3, 1}, {-1, 2}, {-3, -5}, {-7, 1}, {5, 2}, {-7, 0}} {
  }
  virtual ~DerivativesTest() = default;
  // The objects under test
  Gradient<double> tGradDouble;
  Hessian<double> tHessDouble;
  Gradient<std::vector<int>> tGradIntVec;
  Hessian<std::vector<int>> tHessIntVec;
};
/**
 * @test
 * @brief Tests Derivatives.h: Access of elements
 */
TEST_F(DerivativesTest, Access) {
  EXPECT_EQ(tGradDouble.x, -1.0);
  EXPECT_EQ(tGradDouble.y, 2.4);
  EXPECT_EQ(tGradDouble.z, 0.1);
  EXPECT_EQ(tHessDouble.xx, 1.2);
  EXPECT_EQ(tHessDouble.xy, 3.4);
  EXPECT_EQ(tHessDouble.xz, -4.0);
  EXPECT_EQ(tHessDouble.yy, 0.2);
  EXPECT_EQ(tHessDouble.yz, 0.5);
  EXPECT_EQ(tHessDouble.zz, 0.7);
  EXPECT_EQ(tGradIntVec.x[0], 1);
  EXPECT_EQ(tGradIntVec.x[1], 2);
  EXPECT_EQ(tGradIntVec.y[0], -1);
  EXPECT_EQ(tGradIntVec.y[1], -1);
  EXPECT_EQ(tGradIntVec.z[0], 5);
  EXPECT_EQ(tGradIntVec.z[1], 7);
  EXPECT_EQ(tHessIntVec.xx[0], 3);
  EXPECT_EQ(tHessIntVec.xx[1], 1);
  EXPECT_EQ(tHessIntVec.xy[0], -1);
  EXPECT_EQ(tHessIntVec.xy[1], 2);
  EXPECT_EQ(tHessIntVec.xz[0], -3);
  EXPECT_EQ(tHessIntVec.xz[1], -5);
  EXPECT_EQ(tHessIntVec.yy[0], -7);
  EXPECT_EQ(tHessIntVec.yy[1], 1);
  EXPECT_EQ(tHessIntVec.yz[0], 5);
  EXPECT_EQ(tHessIntVec.yz[1], 2);
  EXPECT_EQ(tHessIntVec.zz[0], -7);
  EXPECT_EQ(tHessIntVec.zz[1], 0);
}
/**
 * @test
 * @brief Tests Derivatives.h: Loop over elements with range-based for-loops
 */
TEST_F(DerivativesTest, RangeBasedForLoops) {
  for (auto& component : tGradDouble)
    component *= 2.0;
  EXPECT_EQ(tGradDouble.x, -2.0);
  EXPECT_EQ(tGradDouble.y, 4.8);
  EXPECT_EQ(tGradDouble.z, 0.2);
  for (auto& component : tHessDouble)
    component *= 2.0;
  EXPECT_EQ(tHessDouble.xx, 2.4);
  EXPECT_EQ(tHessDouble.xy, 6.8);
  EXPECT_EQ(tHessDouble.xz, -8.0);
  EXPECT_EQ(tHessDouble.yy, 0.4);
  EXPECT_EQ(tHessDouble.yz, 1.0);
  EXPECT_EQ(tHessDouble.zz, 1.4);
  for (auto& component : tGradIntVec) {
    component[0] *= 2.0;
    component[1] += 1;
  }
  EXPECT_EQ(tGradIntVec.x[0], 2);
  EXPECT_EQ(tGradIntVec.x[1], 3);
  EXPECT_EQ(tGradIntVec.y[0], -2);
  EXPECT_EQ(tGradIntVec.y[1], 0);
  EXPECT_EQ(tGradIntVec.z[0], 10);
  EXPECT_EQ(tGradIntVec.z[1], 8);
  for (auto& component : tHessIntVec) {
    component[0] *= 2.0;
    component[1] += 1;
  }
  EXPECT_EQ(tHessIntVec.xx[0], 6);
  EXPECT_EQ(tHessIntVec.xx[1], 2);
  EXPECT_EQ(tHessIntVec.xy[0], -2);
  EXPECT_EQ(tHessIntVec.xy[1], 3);
  EXPECT_EQ(tHessIntVec.xz[0], -6);
  EXPECT_EQ(tHessIntVec.xz[1], -4);
  EXPECT_EQ(tHessIntVec.yy[0], -14);
  EXPECT_EQ(tHessIntVec.yy[1], 2);
  EXPECT_EQ(tHessIntVec.yz[0], 10);
  EXPECT_EQ(tHessIntVec.yz[1], 3);
  EXPECT_EQ(tHessIntVec.zz[0], -14);
  EXPECT_EQ(tHessIntVec.zz[1], 1);
}
/**
 * @test
 * @brief Tests Derivatives.h: Copy and move construct
 */
TEST_F(DerivativesTest, CopyAndMoveConstruct) {
  Gradient<std::vector<int>> gradCopy(tGradIntVec);
  Hessian<std::vector<int>> hessCopy(tHessIntVec);
  EXPECT_EQ(gradCopy.x[0], 1);
  EXPECT_EQ(gradCopy.x[1], 2);
  EXPECT_EQ(gradCopy.y[0], -1);
  EXPECT_EQ(gradCopy.y[1], -1);
  EXPECT_EQ(gradCopy.z[0], 5);
  EXPECT_EQ(gradCopy.z[1], 7);
  EXPECT_EQ(hessCopy.xx[0], 3);
  EXPECT_EQ(hessCopy.xx[1], 1);
  EXPECT_EQ(hessCopy.xy[0], -1);
  EXPECT_EQ(hessCopy.xy[1], 2);
  EXPECT_EQ(hessCopy.xz[0], -3);
  EXPECT_EQ(hessCopy.xz[1], -5);
  EXPECT_EQ(hessCopy.yy[0], -7);
  EXPECT_EQ(hessCopy.yy[1], 1);
  EXPECT_EQ(hessCopy.yz[0], 5);
  EXPECT_EQ(hessCopy.yz[1], 2);
  EXPECT_EQ(hessCopy.zz[0], -7);
  EXPECT_EQ(hessCopy.zz[1], 0);
  // Move
  Gradient<std::unique_ptr<double>> gradPtr;
  gradPtr.x.reset(new double(5));
  gradPtr.y.reset(new double(4));
  gradPtr.z.reset(new double(1));
  Gradient<std::unique_ptr<double>> gradMove(std::move(gradPtr));
  EXPECT_EQ(*gradMove.x, 5);
  EXPECT_EQ(*gradMove.y, 4);
  EXPECT_EQ(*gradMove.z, 1);
  EXPECT_EQ(gradPtr.x, nullptr);
  EXPECT_EQ(gradPtr.y, nullptr);
  EXPECT_EQ(gradPtr.z, nullptr);
  Hessian<std::unique_ptr<double>> hessPtr;
  hessPtr.xx.reset(new double(1));
  hessPtr.xy.reset(new double(-2.3));
  hessPtr.xz.reset(new double(1.3));
  hessPtr.yy.reset(new double(1.4));
  hessPtr.yz.reset(new double(1.5));
  hessPtr.zz.reset(new double(1.6));
  Hessian<std::unique_ptr<double>> hessMove(std::move(hessPtr));
  EXPECT_EQ(*hessMove.xx, 1);
  EXPECT_EQ(*hessMove.xy, -2.3);
  EXPECT_EQ(*hessMove.xz, 1.3);
  EXPECT_EQ(*hessMove.yy, 1.4);
  EXPECT_EQ(*hessMove.yz, 1.5);
  EXPECT_EQ(*hessMove.zz, 1.6);
  EXPECT_EQ(hessPtr.xx, nullptr);
  EXPECT_EQ(hessPtr.xy, nullptr);
  EXPECT_EQ(hessPtr.xz, nullptr);
  EXPECT_EQ(hessPtr.yy, nullptr);
  EXPECT_EQ(hessPtr.yz, nullptr);
  EXPECT_EQ(hessPtr.zz, nullptr);
}
/**
 * @test
 * @brief Tests Derivatives.h: makeGradient
 */
TEST_F(DerivativesTest, MakeGradient) {
  double foo = 3;
  auto result = makeGradient<double>(foo);
  EXPECT_EQ(result.x, 3);
  EXPECT_EQ(result.y, 3);
  EXPECT_EQ(result.z, 3);
  auto resultPtr = makeGradientPtr<double>(foo);
  EXPECT_EQ(resultPtr->x, 3);
  EXPECT_EQ(resultPtr->y, 3);
  EXPECT_EQ(resultPtr->z, 3);
}
/**
 * @test
 * @brief Tests Derivatives.h: makeHessian
 */
TEST_F(DerivativesTest, MakeHessian) {
  double foo = 3;
  auto result = makeHessian<double>(foo);
  EXPECT_EQ(result.xx, 3);
  EXPECT_EQ(result.xy, 3);
  EXPECT_EQ(result.xz, 3);
  EXPECT_EQ(result.yy, 3);
  EXPECT_EQ(result.yz, 3);
  EXPECT_EQ(result.zz, 3);
  auto resultPtr = makeHessianPtr<double>(foo);
  EXPECT_EQ(resultPtr->xx, 3);
  EXPECT_EQ(resultPtr->xy, 3);
  EXPECT_EQ(resultPtr->xz, 3);
  EXPECT_EQ(resultPtr->yy, 3);
  EXPECT_EQ(resultPtr->yz, 3);
  EXPECT_EQ(resultPtr->zz, 3);
}

} /* namespace Serenity */
