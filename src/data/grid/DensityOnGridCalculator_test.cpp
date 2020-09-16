/**
 * @file   DensityOnGridCalculator_test.cpp
 *
 * @date   Jun 3, 2014
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
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/matrices/DensityMatrix.h"
#include "math/Derivatives.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
#include "testsupply/GridController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class DensityOnGridCalculatorTest : public ::testing::Test {
 protected:
  DensityOnGridCalculatorTest() {
  }
};
/**
 * @test
 * This test relies on a functioning BasisFunctionOnGridController. Reference values were calculated
 * with Mathematica.
 */
TEST_F(DensityOnGridCalculatorTest, TestDensityAndGradientOnGrid) {
  const double expectedPrecision = 5E-8;
  /*
   * Set up test environment
   */
  auto gridController = GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::TINY);
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  /*
   * Set up test environment - Controller
   */
  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController(
      new BasisFunctionOnGridController(basisController, gridController, 3, expectedPrecision * 0.00, 2));
  /*
   * Set up test environment - Test object
   */
  DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED> testObject(basisFunctionOnGridController, 0.0);
  /*
   * Set up test environment - Method input
   */
  DensityMatrix<Options::SCF_MODES::RESTRICTED> densityMatrix(basisController);
  densityMatrix(0, 0) = 1.1 * sqrt(1.1930904 * 1.1930904);
  densityMatrix(0, 1) = -0.15 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(0, 2) = 0.12 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(0, 3) = 0.15 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(0, 4) = 0.02 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(0, 5) = 0.23 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(0, 6) = -0.07 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(0, 7) = 0.23 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(0, 8) = 0.04 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(0, 9) = -0.01 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(1, 1) = 0.25 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(1, 2) = 0.05 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(1, 3) = -0.03 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(1, 4) = 0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(1, 5) = 0.04 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(1, 6) = 0.1 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(1, 7) = 0.07 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(1, 8) = -0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(1, 9) = 0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 2) = 2.31 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(2, 3) = -0.01 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(2, 4) = -0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 5) = 0.1 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 6) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 7) = 0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 8) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 9) = 0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 3) = 0.74 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(3, 4) = 0.12 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 5) = -0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 6) = 0.11 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 7) = 0.003 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 8) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 9) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(4, 4) = 3.13 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(4, 5) = -0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(4, 6) = 0.02 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(4, 7) = -0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(4, 8) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(4, 9) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 5) = 0.66 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 6) = -0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 7) = -0.004 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 8) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 9) = 0.1 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 6) = 0.23 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 7) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 8) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 9) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 7) = 2.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 8) = 0.11 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 9) = -0.005 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 8) = 1.12 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 9) = -0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 9) = 0.11 * sqrt(1.2623580 * 1.2623580);

  densityMatrix(0, 0) = 1.1 * sqrt(1.1930904 * 1.1930904);
  densityMatrix(1, 0) = -0.15 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(2, 0) = 0.12 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(3, 0) = 0.15 * sqrt(1.1930904 * 1.1826193);
  densityMatrix(4, 0) = 0.02 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(5, 0) = 0.23 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(6, 0) = -0.07 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(7, 0) = 0.23 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(8, 0) = 0.04 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(9, 0) = -0.01 * sqrt(1.1930904 * 1.2623580);
  densityMatrix(1, 1) = 0.25 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(2, 1) = 0.05 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(3, 1) = -0.03 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(4, 1) = 0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(5, 1) = 0.04 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(6, 1) = 0.1 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(7, 1) = 0.07 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(8, 1) = -0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(9, 1) = 0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(2, 2) = 2.31 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(3, 2) = -0.01 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(4, 2) = -0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(5, 2) = 0.1 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(6, 2) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(7, 2) = 0.03 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(8, 2) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(9, 2) = 0.02 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(3, 3) = 0.74 * sqrt(1.1826193 * 1.1826193);
  densityMatrix(4, 3) = 0.12 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(5, 3) = -0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(6, 3) = 0.11 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(7, 3) = 0.003 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(8, 3) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(9, 3) = 0.01 * sqrt(1.1826193 * 1.2623580);
  densityMatrix(4, 4) = 3.13 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 4) = -0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 4) = 0.02 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 4) = -0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 4) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 4) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(5, 5) = 0.66 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 5) = -0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 5) = -0.004 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 5) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 5) = 0.1 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(6, 6) = 0.23 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 6) = 0.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 6) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 6) = 0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(7, 7) = 2.01 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 7) = 0.11 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 7) = -0.005 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(8, 8) = 1.12 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 8) = -0.03 * sqrt(1.2623580 * 1.2623580);
  densityMatrix(9, 9) = 0.11 * sqrt(1.2623580 * 1.2623580);
  auto densOnGrid = testObject.calcDensityOnGrid(densityMatrix);
  EXPECT_NEAR(0.998585406945, densOnGrid[0], expectedPrecision);
  EXPECT_NEAR(0.0856147256035, densOnGrid[1], expectedPrecision);
  EXPECT_NEAR(0.159966319953, densOnGrid[2], expectedPrecision);
  EXPECT_NEAR(0.000142599974076, densOnGrid[3], expectedPrecision);
  /*
   * Gradient
   */
  // slightly mess up densOnGrid to make sure it's overwritten
  densOnGrid[0] = 4.32;
  densOnGrid[3] = 1.23490;
  auto gradOnGrid = makeGradient<DensityOnGrid<Options::SCF_MODES::RESTRICTED>>(gridController);
  testObject.calcDensityAndGradientOnGrid(densityMatrix, densOnGrid, gradOnGrid);
  EXPECT_NEAR(0.998585406945, densOnGrid[0], expectedPrecision);
  EXPECT_NEAR(0.0856147256035, densOnGrid[1], expectedPrecision);
  EXPECT_NEAR(0.159966319953, densOnGrid[2], expectedPrecision);
  EXPECT_NEAR(0.000142599974076, densOnGrid[3], expectedPrecision);
  EXPECT_NEAR(0.330967243429, gradOnGrid.x[0], expectedPrecision);
  EXPECT_NEAR(-0.511031515169, gradOnGrid.x[1], expectedPrecision);
  EXPECT_NEAR(0.216024773574, gradOnGrid.x[2], expectedPrecision);
  EXPECT_NEAR(-0.00107982066132, gradOnGrid.x[3], expectedPrecision);
  EXPECT_NEAR(0.109864714691, gradOnGrid.y[0], expectedPrecision);
  EXPECT_NEAR(0.115243549484, gradOnGrid.y[1], expectedPrecision);
  EXPECT_NEAR(0.308780827494, gradOnGrid.y[2], expectedPrecision);
  EXPECT_NEAR(-0.000144049190403, gradOnGrid.y[3], expectedPrecision);
  EXPECT_NEAR(0.173303906447, gradOnGrid.z[0], expectedPrecision);
  EXPECT_NEAR(0.149009688581, gradOnGrid.z[1], expectedPrecision);
  EXPECT_NEAR(-0.0147016342747, gradOnGrid.z[2], expectedPrecision);
  EXPECT_NEAR(0.000480912886676, gradOnGrid.z[3], expectedPrecision);
  /*
   * Hessian
   */
  // slightly mess up old data to make sure it's overwritten
  densOnGrid[1] = -1.02;
  densOnGrid[3] = 99999.9;
  gradOnGrid.x[2] = -12423455;
  gradOnGrid.y[0] = -12423455;
  gradOnGrid.z[1] = -12423455;
  auto hessOnGrid = makeHessian<DensityOnGrid<Options::SCF_MODES::RESTRICTED>>(gridController);
  testObject.calcDensityAndDerivativesOnGrid(densityMatrix, densOnGrid, gradOnGrid, hessOnGrid);
  EXPECT_NEAR(0.998585406945, densOnGrid[0], expectedPrecision);
  EXPECT_NEAR(0.0856147256035, densOnGrid[1], expectedPrecision);
  EXPECT_NEAR(0.159966319953, densOnGrid[2], expectedPrecision);
  EXPECT_NEAR(0.000142599974076, densOnGrid[3], expectedPrecision);
  EXPECT_NEAR(0.330967243429, gradOnGrid.x[0], expectedPrecision);
  EXPECT_NEAR(-0.511031515169, gradOnGrid.x[1], expectedPrecision);
  EXPECT_NEAR(0.216024773574, gradOnGrid.x[2], expectedPrecision);
  EXPECT_NEAR(-0.00107982066132, gradOnGrid.x[3], expectedPrecision);
  EXPECT_NEAR(0.109864714691, gradOnGrid.y[0], expectedPrecision);
  EXPECT_NEAR(0.115243549484, gradOnGrid.y[1], expectedPrecision);
  EXPECT_NEAR(0.308780827494, gradOnGrid.y[2], expectedPrecision);
  EXPECT_NEAR(-0.000144049190403, gradOnGrid.y[3], expectedPrecision);
  EXPECT_NEAR(0.173303906447, gradOnGrid.z[0], expectedPrecision);
  EXPECT_NEAR(0.149009688581, gradOnGrid.z[1], expectedPrecision);
  EXPECT_NEAR(-0.0147016342747, gradOnGrid.z[2], expectedPrecision);
  EXPECT_NEAR(0.000480912886676, gradOnGrid.z[3], expectedPrecision);
  EXPECT_NEAR(-4.0497395077622, hessOnGrid.xx[0], expectedPrecision);
  EXPECT_NEAR(3.3946178224277, hessOnGrid.xx[1], expectedPrecision);
  EXPECT_NEAR(0.17663869285660, hessOnGrid.xx[2], expectedPrecision);
  EXPECT_NEAR(0.0076116193299364, hessOnGrid.xx[3], expectedPrecision);
  EXPECT_NEAR(0.077547198866445, hessOnGrid.xy[0], expectedPrecision);
  EXPECT_NEAR(0.058695288791396, hessOnGrid.xy[1], expectedPrecision);
  EXPECT_NEAR(0.30516290080961, hessOnGrid.xy[2], expectedPrecision);
  EXPECT_NEAR(0.0010355428268474, hessOnGrid.xy[3], expectedPrecision);
  EXPECT_NEAR(0.34553903256919, hessOnGrid.xz[0], expectedPrecision);
  EXPECT_NEAR(-0.50017522084132, hessOnGrid.xz[1], expectedPrecision);
  EXPECT_NEAR(-0.035922086723103, hessOnGrid.xz[2], expectedPrecision);
  EXPECT_NEAR(-0.0036084520539767, hessOnGrid.xz[3], expectedPrecision);
  EXPECT_NEAR(-2.9308731686653, hessOnGrid.yy[0], expectedPrecision);
  EXPECT_NEAR(14.009221944124, hessOnGrid.yy[1], expectedPrecision);
  EXPECT_NEAR(0.38202375529246, hessOnGrid.yy[2], expectedPrecision);
  EXPECT_NEAR(-0.00030798024947144, hessOnGrid.yy[3], expectedPrecision);
  EXPECT_NEAR(0.022894432435487, hessOnGrid.yz[0], expectedPrecision);
  EXPECT_NEAR(-0.061730533541120, hessOnGrid.yz[1], expectedPrecision);
  EXPECT_NEAR(-0.026389142078342, hessOnGrid.yz[2], expectedPrecision);
  EXPECT_NEAR(-0.00039799531131632, hessOnGrid.yz[3], expectedPrecision);
  EXPECT_NEAR(-4.0126553586185, hessOnGrid.zz[0], expectedPrecision);
  EXPECT_NEAR(4.2523888726439, hessOnGrid.zz[1], expectedPrecision);
  EXPECT_NEAR(-0.23944457823680, hessOnGrid.zz[2], expectedPrecision);
  EXPECT_NEAR(0.0010780480705487, hessOnGrid.zz[3], expectedPrecision);
}

} // namespace Serenity
