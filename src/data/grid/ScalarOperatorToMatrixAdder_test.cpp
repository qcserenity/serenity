/**
 * @file   ScalarOperatorToMatrixAdder_test.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   23. September 2015, 16:09
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
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/GridPotential.h"
#include "data/matrices/FockMatrix.h"
#include "math/Derivatives.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
#include "testsupply/GridController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ScalarOperatorToMatrixAdderTest : public testing::Test {
 protected:
};
/**
 * @test
 * This test relies on a functioning BasisFunctionOnGridController. Reference values were calculated
 * with Mathematica.
 */
TEST_F(ScalarOperatorToMatrixAdderTest, IntegrationOfScalarAndGradientPotentialsRestricted) {
  /*
   * TODO original accuracy was 1E-12, but to respect a change in the normalization of the basis
   * functions the results were modified by less accurate numbers. A recalculation of the values
   * with Mathematica is needed.
   */
  const double expectedPrecision = 1E-8;
  /*
   * Set up test environment
   */
  auto gridController = GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::VERY_SMALL);
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  /*
   * Set up test environment - Controller
   */
  std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController(
      new BasisFunctionOnGridController(basisController, gridController, 3, expectedPrecision * 1e-3, 1));
  /*
   * Test object
   */
  ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED> rTestObject(basisFunctionOnGridController,
                                                                          expectedPrecision * 1e-3);
  /*
   * Define input data
   */
  GridPotential<Options::SCF_MODES::RESTRICTED> pot(gridController);
  pot[0] = 1.3;
  pot[1] = -0.7;
  pot[2] = 5.4;
  pot[3] = 0.0003;
  pot[4] = -0.3;
  auto gradPot = makeGradient<GridPotential<Options::SCF_MODES::RESTRICTED>>(gridController);
  gradPot.x[0] = -1.0 * -0.0002;
  gradPot.x[1] = -1.0 * 0.02;
  gradPot.x[2] = -1.0 * -1.3;
  gradPot.x[3] = -1.0 * -0.004;
  gradPot.x[4] = -1.0 * 0.001;
  gradPot.y[0] = -1.0 * 0.01;
  gradPot.y[1] = -1.0 * -0.1;
  gradPot.y[2] = -1.0 * 2.2;
  gradPot.y[3] = -1.0 * 0.03;
  gradPot.y[4] = -1.0 * -0.0025;
  gradPot.z[0] = -1.0 * -1.03;
  gradPot.z[1] = -1.0 * 0.007;
  gradPot.z[2] = -1.0 * -0.4;
  gradPot.z[3] = -1.0 * 0.0000007;
  gradPot.z[4] = -1.0 * 0.00004;
  /*
   * Calculate and check result
   */
  FockMatrix<Options::SCF_MODES::RESTRICTED> result(basisController);
  rTestObject.addScalarOperatorToMatrix(result, pot, gradPot);
  EXPECT_NEAR(result(0, 0), -0.014723532368381934, expectedPrecision);
  EXPECT_NEAR(result(0, 1), -0.016771151488521534, expectedPrecision);
  EXPECT_NEAR(result(0, 2), -0.1132735667955283, expectedPrecision);
  EXPECT_NEAR(result(0, 3), -0.071518029901160199, expectedPrecision);
  EXPECT_NEAR(result(0, 4), -0.0020107774569296821, expectedPrecision);
  EXPECT_NEAR(result(0, 5), 0.0077303812833758659, expectedPrecision);
  EXPECT_NEAR(result(0, 6), 0.0025490282389903883, expectedPrecision);
  EXPECT_NEAR(result(0, 7), -0.0066319996220955335, expectedPrecision);
  EXPECT_NEAR(result(0, 8), -0.0047794145660672279, expectedPrecision);
  EXPECT_NEAR(result(0, 9), -0.0013473580158660139, expectedPrecision);
  EXPECT_NEAR(result(1, 1), 0.032751588205558152, expectedPrecision);
  EXPECT_NEAR(result(1, 2), 0.090170891284227783, expectedPrecision);
  EXPECT_NEAR(result(1, 3), 0.064062754837144284, expectedPrecision);
  EXPECT_NEAR(result(1, 4), 0.0014053014611159543, expectedPrecision);
  EXPECT_NEAR(result(1, 5), -0.0052954871245275816, expectedPrecision);
  EXPECT_NEAR(result(1, 6), -0.0016737378844947513, expectedPrecision);
  EXPECT_NEAR(result(1, 7), 0.0061180244213455623, expectedPrecision);
  EXPECT_NEAR(result(1, 8), 0.0036458070154504153, expectedPrecision);
  EXPECT_NEAR(result(1, 9), 0.00083920101159597887, expectedPrecision);
  EXPECT_NEAR(result(2, 2), 0.046423374970203837, expectedPrecision);
  EXPECT_NEAR(result(2, 3), 0.041548440855128294, expectedPrecision);
  EXPECT_NEAR(result(2, 4), -0.0005591466884157651, expectedPrecision);
  EXPECT_NEAR(result(2, 5), 0.0021233203914992099, expectedPrecision);
  EXPECT_NEAR(result(2, 6), 0.00070437755445713196, expectedPrecision);
  EXPECT_NEAR(result(2, 7), -0.0027622566895728552, expectedPrecision);
  EXPECT_NEAR(result(2, 8), -0.0016130435192994251, expectedPrecision);
  EXPECT_NEAR(result(2, 9), -0.00032710416667275973, expectedPrecision);
  EXPECT_NEAR(result(3, 3), 0.023426478769345564, expectedPrecision);
  EXPECT_NEAR(result(3, 4), -0.00052168170493647613, expectedPrecision);
  EXPECT_NEAR(result(3, 5), 0.001905187547403903, expectedPrecision);
  EXPECT_NEAR(result(3, 6), 0.0004855463346554416, expectedPrecision);
  EXPECT_NEAR(result(3, 7), -0.0021211696229092446, expectedPrecision);
  EXPECT_NEAR(result(3, 8), -0.00082728519755134982, expectedPrecision);
  EXPECT_NEAR(result(3, 9), -4.6931474442151502e-05, expectedPrecision);
  EXPECT_NEAR(result(4, 4), -4.5055793961758299e-05, expectedPrecision);
  EXPECT_NEAR(result(4, 5), 0.00017118796018732356, expectedPrecision);
  EXPECT_NEAR(result(4, 6), 4.5390973760116205e-05, expectedPrecision);
  EXPECT_NEAR(result(4, 7), -0.00021703940343490759, expectedPrecision);
  EXPECT_NEAR(result(4, 8), -0.00010200506291916221, expectedPrecision);
  EXPECT_NEAR(result(4, 9), -1.6497218574881483e-05, expectedPrecision);
  EXPECT_NEAR(result(5, 5), -0.00065111821030472265, expectedPrecision);
  EXPECT_NEAR(result(5, 6), -0.00017667795160524906, expectedPrecision);
  EXPECT_NEAR(result(5, 7), 0.00082679829543728077, expectedPrecision);
  EXPECT_NEAR(result(5, 8), 0.00039929366100157682, expectedPrecision);
  EXPECT_NEAR(result(5, 9), 6.7355953660195146e-05, expectedPrecision);
  EXPECT_NEAR(result(6, 6), -4.9491655724644462e-05, expectedPrecision);
  EXPECT_NEAR(result(6, 7), 0.00023053230266497154, expectedPrecision);
  EXPECT_NEAR(result(6, 8), 0.00011666393393171289, expectedPrecision);
  EXPECT_NEAR(result(6, 9), 2.1586419959790338e-05, expectedPrecision);
  EXPECT_NEAR(result(7, 7), 0.00011324435021710181, expectedPrecision);
  EXPECT_NEAR(result(7, 8), 3.4779152585474543e-05, expectedPrecision);
  EXPECT_NEAR(result(7, 9), -2.1029250518470563e-05, expectedPrecision);
  EXPECT_NEAR(result(8, 8), -6.3087751555411655e-05, expectedPrecision);
  EXPECT_NEAR(result(8, 9), -4.493481649631375e-05, expectedPrecision);
  EXPECT_NEAR(result(9, 9), -1.936970873014834e-05, expectedPrecision);
};

} /* namespace Serenity */
