/**
 * @file   BasisFunctionOnGridController_test.cpp
 *
 * @date   Jun 3, 2014, last rework June 6, 2017
 * @author Thomas Dresselhaus, last rework Jan Unsleber
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
#include "data/grid/BasisFunctionOnGridController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/Normalization.h"
#include "math/Derivatives.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
#include "testsupply/GridController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class BasisFunctionOnGridControllerTest : public testing::Test {};
/**
 * @test
 * @brief Tested against values which were calculated manually.
 *
 * Results were calculated with Mathematica. Its analytic derivatives were used.
 */
TEST_F(BasisFunctionOnGridControllerTest, TestBasisFunctionValuesAndDerivatives) {
  /*
   * TODO original accuracy was 1E-11, but to respect a change in the normalization of the basis
   * functions the results were modified by less accurate numbers. A recalculation of the values
   * with Mathematica is needed.
   */
  const double expectedPrecision = 1E-8;
  /*
   * Set up test environment
   */
  auto gridController = GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::TINY);
  auto basisController = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  /*
   * Set up test environment - Test object
   */
  BasisFunctionOnGridController testObject(basisController, gridController, 10, expectedPrecision * 0.001, 2);

  /*
   * Test methods
   */
  EXPECT_EQ((unsigned int)1, testObject.getNBlocks());
  EXPECT_EQ((unsigned int)4, testObject.getNGridPoints());
  EXPECT_EQ((unsigned int)0, testObject.getFirstIndexOfBlock(0));
  auto& basFuncData = testObject.getBlockOnGridData(0);
  const Eigen::MatrixXd& values = basFuncData->functionValues;
  ASSERT_EQ(values.cols(), basisController->getNBasisFunctions());
  ASSERT_EQ(values.rows(), (unsigned int)4);
  const auto& derivatives = *basFuncData->derivativeValues;
  for (auto& component : derivatives) {
    ASSERT_EQ(component.cols(), basisController->getNBasisFunctions());
    ASSERT_EQ(component.rows(), (unsigned int)4);
  }
  const auto& hessian = *basFuncData->secondDerivativeValues;
  for (auto& component : hessian) {
    ASSERT_EQ(component.cols(), basisController->getNBasisFunctions());
    ASSERT_EQ(component.rows(), (unsigned int)4);
  }

  /*
   * Check the calculated basis function values against manually calculated references
   */
  EXPECT_NEAR(0.832567765386 / sqrt(1.1930904), values(0, 0), expectedPrecision);
  EXPECT_NEAR(0.278411287801 / sqrt(1.1930904), values(1, 0), expectedPrecision);
  EXPECT_NEAR(0.00176735669536 / sqrt(1.1930904), values(2, 0), expectedPrecision);
  EXPECT_NEAR(0.0000251996341204 / sqrt(1.1930904), values(3, 0), expectedPrecision);
  EXPECT_NEAR(-0.570260987101 / sqrt(1.1826193), values(0, 1), expectedPrecision);
  EXPECT_NEAR(0, values(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.0000714206448118 / sqrt(1.1826193), values(2, 1), expectedPrecision);
  EXPECT_NEAR(0.0149784111869 / sqrt(1.1826193), values(3, 1), expectedPrecision);
  EXPECT_NEAR(0, values(0, 2), expectedPrecision);
  EXPECT_NEAR(0, values(1, 2), expectedPrecision);
  EXPECT_NEAR(0.0000238068816039 / sqrt(1.1826193), values(2, 2), expectedPrecision);
  EXPECT_NEAR(0.00374460279674 / sqrt(1.1826193), values(3, 2), expectedPrecision);
  EXPECT_NEAR(0, values(0, 3), expectedPrecision);
  EXPECT_NEAR(0, values(1, 3), expectedPrecision);
  EXPECT_NEAR(0.0000238068816039 / sqrt(1.1826193), values(2, 3), expectedPrecision);
  EXPECT_NEAR(-0.00748920559347 / sqrt(1.1826193), values(3, 3), expectedPrecision);
  EXPECT_NEAR(0.00637294953921 / sqrt(1.2623580), values(0, 4), expectedPrecision);
  EXPECT_NEAR(0.00565239243008 / sqrt(1.2623580), values(1, 4), expectedPrecision);
  EXPECT_NEAR(0.0954200464484 / sqrt(1.2623580), values(2, 4), expectedPrecision);
  EXPECT_NEAR(0.0000384842654340 / sqrt(1.2623580), values(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.0275956809900 / sqrt(1.2623580), values(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0122377885915 / sqrt(1.2623580), values(1, 5), expectedPrecision);
  EXPECT_NEAR(0.247908552764 / sqrt(1.2623580), values(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000333283515118 / sqrt(1.2623580), values(3, 5), expectedPrecision);
  EXPECT_NEAR(-0.0110382723960 / sqrt(1.2623580), values(0, 6), expectedPrecision);
  EXPECT_NEAR(-0.00489511543661 / sqrt(1.2623580), values(1, 6), expectedPrecision);
  EXPECT_NEAR(0, values(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000333283515118 / sqrt(1.2623580), values(3, 6), expectedPrecision);
  EXPECT_NEAR(0.0398309346201 / sqrt(1.2623580), values(0, 7), expectedPrecision);
  EXPECT_NEAR(0.00883186317200 / sqrt(1.2623580), values(1, 7), expectedPrecision);
  EXPECT_NEAR(0.214695104509 / sqrt(1.2623580), values(2, 7), expectedPrecision);
  EXPECT_NEAR(9.62106635849e-6 / sqrt(1.2623580), values(3, 7), expectedPrecision);
  EXPECT_NEAR(0.0275956809900 / sqrt(1.2623580), values(0, 8), expectedPrecision);
  EXPECT_NEAR(0.00611889429576 / sqrt(1.2623580), values(1, 8), expectedPrecision);
  EXPECT_NEAR(0, values(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000166641757559 / sqrt(1.2623580), values(3, 8), expectedPrecision);
  EXPECT_NEAR(0.00637294953921 / sqrt(1.2623580), values(0, 9), expectedPrecision);
  EXPECT_NEAR(0.00141309810752 / sqrt(1.2623580), values(1, 9), expectedPrecision);
  EXPECT_NEAR(0, values(2, 9), expectedPrecision);
  EXPECT_NEAR(9.62106635849e-6 / sqrt(1.2623580), values(3, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(0, 0), expectedPrecision);
  EXPECT_NEAR(-0.589265770897 / sqrt(1.1930904), derivatives.x(1, 0), expectedPrecision);
  EXPECT_NEAR(0.00707237261901 / sqrt(1.1930904), derivatives.x(2, 0), expectedPrecision);
  EXPECT_NEAR(-0.000151198703799 / sqrt(1.1930904), derivatives.x(3, 0), expectedPrecision);
  EXPECT_NEAR(-0.662024200686 / sqrt(1.1826193), derivatives.x(0, 1), expectedPrecision);
  EXPECT_NEAR(1.76443270721 / sqrt(1.1826193), derivatives.x(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.000404718689509 / sqrt(1.1826193), derivatives.x(2, 1), expectedPrecision);
  EXPECT_NEAR(-0.0524991227999 / sqrt(1.1826193), derivatives.x(3, 1), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(0, 2), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(1, 2), expectedPrecision);
  EXPECT_NEAR(0.000142841857038 / sqrt(1.1826193), derivatives.x(2, 2), expectedPrecision);
  EXPECT_NEAR(-0.0149970820983 / sqrt(1.1826193), derivatives.x(3, 2), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(0, 3), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(1, 3), expectedPrecision);
  EXPECT_NEAR(0.000142841857038 / sqrt(1.1826193), derivatives.x(2, 3), expectedPrecision);
  EXPECT_NEAR(0.0299941641967 / sqrt(1.1826193), derivatives.x(3, 3), expectedPrecision);
  EXPECT_NEAR(0.00630932138969 / sqrt(1.2623580), derivatives.x(0, 4), expectedPrecision);
  EXPECT_NEAR(-0.00566630089326 / sqrt(1.2623580), derivatives.x(1, 4), expectedPrecision);
  EXPECT_NEAR(-0.0697506188560 / sqrt(1.2623580), derivatives.x(2, 4), expectedPrecision);
  EXPECT_NEAR(-0.000134694954234 / sqrt(1.2623580), derivatives.x(3, 4), expectedPrecision);
  EXPECT_NEAR(0.000275517969394 / sqrt(1.2623580), derivatives.x(0, 5), expectedPrecision);
  EXPECT_NEAR(0.0183867955934 / sqrt(1.2623580), derivatives.x(1, 5), expectedPrecision);
  EXPECT_NEAR(0.0666911291870 / sqrt(1.2623580), derivatives.x(2, 5), expectedPrecision);
  EXPECT_NEAR(0.000124981340006 / sqrt(1.2623580), derivatives.x(3, 5), expectedPrecision);
  EXPECT_NEAR(0.000110207187758 / sqrt(1.2623580), derivatives.x(0, 6), expectedPrecision);
  EXPECT_NEAR(0.00735471823735 / sqrt(1.2623580), derivatives.x(1, 6), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(2, 6), expectedPrecision);
  EXPECT_NEAR(0.000124981340006 / sqrt(1.2623580), derivatives.x(3, 6), expectedPrecision);
  EXPECT_NEAR(-0.0402286105545 / sqrt(1.2623580), derivatives.x(0, 7), expectedPrecision);
  EXPECT_NEAR(-0.0176854583177 / sqrt(1.2623580), derivatives.x(1, 7), expectedPrecision);
  EXPECT_NEAR(0.272451316592 / sqrt(1.2623580), derivatives.x(2, 7), expectedPrecision);
  EXPECT_NEAR(-0.0000384842717378 / sqrt(1.2623580), derivatives.x(3, 7), expectedPrecision);
  EXPECT_NEAR(-0.0278711989593 / sqrt(1.2623580), derivatives.x(0, 8), expectedPrecision);
  EXPECT_NEAR(-0.0122528449446 / sqrt(1.2623580), derivatives.x(1, 8), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(2, 8), expectedPrecision);
  EXPECT_NEAR(-0.0000666567139421 / sqrt(1.2623580), derivatives.x(3, 8), expectedPrecision);
  EXPECT_NEAR(-0.00643657768873 / sqrt(1.2623580), derivatives.x(0, 9), expectedPrecision);
  EXPECT_NEAR(-0.00282967333083 / sqrt(1.2623580), derivatives.x(1, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.x(2, 9) / sqrt(1.2623580), expectedPrecision);
  EXPECT_NEAR(-0.0000384842717378 / sqrt(1.2623580), derivatives.x(3, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(0, 0), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(1, 0), expectedPrecision);
  EXPECT_NEAR(-0.00353618630950 / sqrt(1.1930904), derivatives.y(2, 0), expectedPrecision);
  EXPECT_NEAR(-0.0000251997839666 / sqrt(1.1930904), derivatives.y(3, 0), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(0, 1), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(1, 1), expectedPrecision);
  EXPECT_NEAR(0.000142841857038 / sqrt(1.1826193), derivatives.y(2, 1), expectedPrecision);
  EXPECT_NEAR(-0.0149970820983 / sqrt(1.1826193), derivatives.y(3, 1), expectedPrecision);
  EXPECT_NEAR(0.570260987101 / sqrt(1.1826193), derivatives.y(0, 2), expectedPrecision);
  EXPECT_NEAR(1.76443270721 / sqrt(1.1826193), derivatives.y(1, 2), expectedPrecision);
  EXPECT_NEAR(-0.0000238070707420 / sqrt(1.1826193), derivatives.y(2, 2), expectedPrecision);
  EXPECT_NEAR(0.00373993506889 / sqrt(1.1826193), derivatives.y(3, 2), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(0, 3), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(1, 3), expectedPrecision);
  EXPECT_NEAR(-0.0000476139523459 / sqrt(1.1826193), derivatives.y(2, 3), expectedPrecision);
  EXPECT_NEAR(0.00749854104917 / sqrt(1.1826193), derivatives.y(3, 3), expectedPrecision);
  EXPECT_NEAR(0.0160914442218 / sqrt(1.2623580), derivatives.y(0, 4), expectedPrecision);
  EXPECT_NEAR(0.0141483666542 / sqrt(1.2623580), derivatives.y(1, 4), expectedPrecision);
  EXPECT_NEAR(0.181634211061 / sqrt(1.2623580), derivatives.y(2, 4), expectedPrecision);
  EXPECT_NEAR(0.0000769685434756 / sqrt(1.2623580), derivatives.y(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.0586397250024 / sqrt(1.2623580), derivatives.y(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0257369969248 / sqrt(1.2623580), derivatives.y(1, 5), expectedPrecision);
  EXPECT_NEAR(0.306627154417 / sqrt(1.2623580), derivatives.y(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000499925381862 / sqrt(1.2623580), derivatives.y(3, 5), expectedPrecision);
  EXPECT_NEAR(-0.0278711989593 / sqrt(1.2623580), derivatives.y(0, 6), expectedPrecision);
  EXPECT_NEAR(-0.0122528449446 / sqrt(1.2623580), derivatives.y(1, 6), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000666567139421 / sqrt(1.2623580), derivatives.y(3, 6), expectedPrecision);
  EXPECT_NEAR(0.0687067786903 / sqrt(1.2623580), derivatives.y(0, 7), expectedPrecision);
  EXPECT_NEAR(0.0150413323595 / sqrt(1.2623580), derivatives.y(1, 7), expectedPrecision);
  EXPECT_NEAR(0.122416835543 / sqrt(1.2623580), derivatives.y(2, 7), expectedPrecision);
  EXPECT_NEAR(9.62106951040e-6 / sqrt(1.2623580), derivatives.y(3, 7), expectedPrecision);
  EXPECT_NEAR(0.0586397250024 / sqrt(1.2623580), derivatives.y(0, 8), expectedPrecision);
  EXPECT_NEAR(0.0128684984624 / sqrt(1.2623580), derivatives.y(1, 8), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000249962690931 / sqrt(1.2623580), derivatives.y(3, 8), expectedPrecision);
  EXPECT_NEAR(0.0160914442218 / sqrt(1.2623580), derivatives.y(0, 9), expectedPrecision);
  EXPECT_NEAR(0.00353709166354 / sqrt(1.2623580), derivatives.y(1, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.y(2, 9), expectedPrecision);
  EXPECT_NEAR(0.0000192421358689 / sqrt(1.2623580), derivatives.y(3, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(0, 0), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(1, 0), expectedPrecision);
  EXPECT_NEAR(-0.00353618630950 / sqrt(1.1930904), derivatives.z(2, 0), expectedPrecision);
  EXPECT_NEAR(0.0000503995679331 / sqrt(1.1930904), derivatives.z(3, 0), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(0, 1), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(1, 1), expectedPrecision);
  EXPECT_NEAR(0.000142841857038 / sqrt(1.1826193), derivatives.z(2, 1), expectedPrecision);
  EXPECT_NEAR(0.0299941641967 / sqrt(1.1826193), derivatives.z(3, 1), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(0, 2), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(1, 2), expectedPrecision);
  EXPECT_NEAR(-0.0000476139523459 / sqrt(1.1826193), derivatives.z(2, 2), expectedPrecision);
  EXPECT_NEAR(0.00749854104917 / sqrt(1.1826193), derivatives.z(3, 2), expectedPrecision);
  EXPECT_NEAR(0.570260987101 / sqrt(1.1826193), derivatives.z(0, 3), expectedPrecision);
  EXPECT_NEAR(1.76443270721 / sqrt(1.1826193), derivatives.z(1, 3), expectedPrecision);
  EXPECT_NEAR(-0.0000238070707420 / sqrt(1.1826193), derivatives.z(2, 3), expectedPrecision);
  EXPECT_NEAR(-0.00750787650488 / sqrt(1.1826193), derivatives.z(3, 3), expectedPrecision);
  EXPECT_NEAR(0.00643657768873 / sqrt(1.2623580), derivatives.z(0, 4), expectedPrecision);
  EXPECT_NEAR(0.00565934666167 / sqrt(1.2623580), derivatives.z(1, 4), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(2, 4), expectedPrecision);
  EXPECT_NEAR(0.0000769685434756 / sqrt(1.2623580), derivatives.z(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.0278711989593 / sqrt(1.2623580), derivatives.z(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0122528449446 / sqrt(1.2623580), derivatives.z(1, 5), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000666567139421 / sqrt(1.2623580), derivatives.z(3, 5), expectedPrecision);
  EXPECT_NEAR(-0.000110207187758 / sqrt(1.2623580), derivatives.z(0, 6), expectedPrecision);
  EXPECT_NEAR(-6.02254121998e-6 / sqrt(1.2623580), derivatives.z(1, 6), expectedPrecision);
  EXPECT_NEAR(-0.165272368509 / sqrt(1.2623580), derivatives.z(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000499925381862 / sqrt(1.2623580), derivatives.z(3, 6), expectedPrecision);
  EXPECT_NEAR(0.0402286105545 / sqrt(1.2623580), derivatives.z(0, 7), expectedPrecision);
  EXPECT_NEAR(0.00884272915885 / sqrt(1.2623580), derivatives.z(1, 7), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(2, 7), expectedPrecision);
  EXPECT_NEAR(0.0000192421358689 / sqrt(1.2623580), derivatives.z(3, 7), expectedPrecision);
  EXPECT_NEAR(0.000275517969394 / sqrt(1.2623580), derivatives.z(0, 8), expectedPrecision);
  EXPECT_NEAR(7.52817652497e-6 / sqrt(1.2623580), derivatives.z(1, 8), expectedPrecision);
  EXPECT_NEAR(-0.247908552764 / sqrt(1.2623580), derivatives.z(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000249962690931 / sqrt(1.2623580), derivatives.z(3, 8), expectedPrecision);
  EXPECT_NEAR(-0.00630932138969 / sqrt(1.2623580), derivatives.z(0, 9), expectedPrecision);
  EXPECT_NEAR(-0.00141135954962 / sqrt(1.2623580), derivatives.z(1, 9), expectedPrecision);
  EXPECT_NEAR(0, derivatives.z(2, 9), expectedPrecision);
  EXPECT_NEAR(9.62106951040e-6 / sqrt(1.2623580), derivatives.z(3, 9), expectedPrecision);
  EXPECT_NEAR(-1.90486012083 / sqrt(1.1930904), hessian.xx(0, 0), expectedPrecision);
  EXPECT_NEAR(0.719038552076 / sqrt(1.1930904), hessian.xx(1, 0), expectedPrecision);
  EXPECT_NEAR(0.0247768708671 / sqrt(1.1930904), hessian.xx(2, 0), expectedPrecision);
  EXPECT_NEAR(0.000856803443789 / sqrt(1.1930904), hessian.xx(3, 0), expectedPrecision);
  EXPECT_NEAR(0.865232333446 / sqrt(1.1826193), hessian.xx(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.00214264828248 / sqrt(1.1826193), hessian.xx(2, 1), expectedPrecision);
  EXPECT_NEAR(0.150568290148 / sqrt(1.1826193), hessian.xx(3, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(1, 2), expectedPrecision);
  EXPECT_NEAR(0.000809443998852 / sqrt(1.1826193), hessian.xx(2, 2), expectedPrecision);
  EXPECT_NEAR(0.0526391546354 / sqrt(1.1826193), hessian.xx(3, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(1, 3), expectedPrecision);
  EXPECT_NEAR(0.000809443998852 / sqrt(1.1826193), hessian.xx(2, 3), expectedPrecision);
  EXPECT_NEAR(-0.105278309271 / sqrt(1.1826193), hessian.xx(3, 3), expectedPrecision);
  EXPECT_NEAR(-0.0128477041176 / sqrt(1.2623580), hessian.xx(0, 4), expectedPrecision);
  EXPECT_NEAR(-0.00276638982337 / sqrt(1.2623580), hessian.xx(1, 4), expectedPrecision);
  EXPECT_NEAR(-0.231911177045 / sqrt(1.2623580), hessian.xx(2, 4), expectedPrecision);
  EXPECT_NEAR(0.000428137764361 / sqrt(1.2623580), hessian.xx(3, 4), expectedPrecision);
  EXPECT_NEAR(0.0550811547921 / sqrt(1.2623580), hessian.xx(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0123973859338 / sqrt(1.2623580), hessian.xx(1, 5), expectedPrecision);
  EXPECT_NEAR(-0.469140653853 / sqrt(1.2623580), hessian.xx(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.000433268850260 / sqrt(1.2623580), hessian.xx(3, 5), expectedPrecision);
  EXPECT_NEAR(0.0220324619169 / sqrt(1.2623580), hessian.xx(0, 6), expectedPrecision);
  EXPECT_NEAR(-0.00495895437354 / sqrt(1.2623580), hessian.xx(1, 6), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.000433268850260 / sqrt(1.2623580), hessian.xx(3, 6), expectedPrecision);
  EXPECT_NEAR(0.000954422242778 / sqrt(1.2623580), hessian.xx(0, 7), expectedPrecision);
  EXPECT_NEAR(0.0266325009504 / sqrt(1.2623580), hessian.xx(1, 7), expectedPrecision);
  EXPECT_NEAR(0.138614908999 / sqrt(1.2623580), hessian.xx(2, 7), expectedPrecision);
  EXPECT_NEAR(0.000144316079533 / sqrt(1.2623580), hessian.xx(3, 7), expectedPrecision);
  EXPECT_NEAR(0.000661243126546 / sqrt(1.2623580), hessian.xx(0, 8), expectedPrecision);
  EXPECT_NEAR(0.0184515379115 / sqrt(1.2623580), hessian.xx(1, 8), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(2, 8), expectedPrecision);
  EXPECT_NEAR(0.000249962782101 / sqrt(1.2623580), hessian.xx(3, 8), expectedPrecision);
  EXPECT_NEAR(0.000152707558845 / sqrt(1.2623580), hessian.xx(0, 9), expectedPrecision);
  EXPECT_NEAR(0.00426120015206 / sqrt(1.2623580), hessian.xx(1, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.xx(2, 9), expectedPrecision);
  EXPECT_NEAR(0.000144316079533 / sqrt(1.2623580), hessian.xx(3, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(0, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(1, 0), expectedPrecision);
  EXPECT_NEAR(-0.0141565285883 / sqrt(1.1930904), hessian.xy(2, 0), expectedPrecision);
  EXPECT_NEAR(0.000151200501954 / sqrt(1.1930904), hessian.xy(3, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(1, 1), expectedPrecision);
  EXPECT_NEAR(0.000809443998852 / sqrt(1.1826193), hessian.xy(2, 1), expectedPrecision);
  EXPECT_NEAR(0.0526391546354 / sqrt(1.1826193), hessian.xy(3, 1), expectedPrecision);
  EXPECT_NEAR(1.23228518779 / sqrt(1.1826193), hessian.xy(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(1, 2), expectedPrecision);
  EXPECT_NEAR(-0.000142844126695 / sqrt(1.1826193), hessian.xy(2, 2), expectedPrecision);
  EXPECT_NEAR(-0.0149597402755 / sqrt(1.1826193), hessian.xy(3, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(1, 3), expectedPrecision);
  EXPECT_NEAR(-0.000285685983732 / sqrt(1.1826193), hessian.xy(2, 3), expectedPrecision);
  EXPECT_NEAR(-0.0300688478423 / sqrt(1.1826193), hessian.xy(3, 3), expectedPrecision);
  EXPECT_NEAR(0.0157096753247 / sqrt(1.2623580), hessian.xy(0, 4), expectedPrecision);
  EXPECT_NEAR(-0.0142318174332 / sqrt(1.2623580), hessian.xy(1, 4), expectedPrecision);
  EXPECT_NEAR(-0.0892242717285 / sqrt(1.2623580), hessian.xy(2, 4), expectedPrecision);
  EXPECT_NEAR(-0.000269390023198 / sqrt(1.2623580), hessian.xy(3, 4), expectedPrecision);
  EXPECT_NEAR(0.00154290062861 / sqrt(1.2623580), hessian.xy(0, 5), expectedPrecision);
  EXPECT_NEAR(0.0387741265414 / sqrt(1.2623580), hessian.xy(1, 5), expectedPrecision);
  EXPECT_NEAR(0.195627312282 / sqrt(1.2623580), hessian.xy(2, 5), expectedPrecision);
  EXPECT_NEAR(0.000187472112098 / sqrt(1.2623580), hessian.xy(3, 5), expectedPrecision);
  EXPECT_NEAR(0.000661243126546 / sqrt(1.2623580), hessian.xy(0, 6), expectedPrecision);
  EXPECT_NEAR(0.0184515379115 / sqrt(1.2623580), hessian.xy(1, 6), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(2, 6), expectedPrecision);
  EXPECT_NEAR(0.000249962782101 / sqrt(1.2623580), hessian.xy(3, 6), expectedPrecision);
  EXPECT_NEAR(-0.0707746935497 / sqrt(1.2623580), hessian.xy(0, 7), expectedPrecision);
  EXPECT_NEAR(-0.0301956709824 / sqrt(1.2623580), hessian.xy(1, 7), expectedPrecision);
  EXPECT_NEAR(0.253330916264 / sqrt(1.2623580), hessian.xy(2, 7), expectedPrecision);
  EXPECT_NEAR(-0.0000384843019961 / sqrt(1.2623580), hessian.xy(3, 7), expectedPrecision);
  EXPECT_NEAR(-0.0601826256310 / sqrt(1.2623580), hessian.xy(0, 8), expectedPrecision);
  EXPECT_NEAR(-0.0258213125019 / sqrt(1.2623580), hessian.xy(1, 8), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(2, 8), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.xy(3, 8), expectedPrecision);
  EXPECT_NEAR(-0.0164732131189 / sqrt(1.2623580), hessian.xy(0, 9), expectedPrecision);
  EXPECT_NEAR(-0.00709504602185 / sqrt(1.2623580), hessian.xy(1, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.xy(2, 9), expectedPrecision);
  EXPECT_NEAR(-0.0000769685737339 / sqrt(1.2623580), hessian.xy(3, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(0, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(1, 0), expectedPrecision);
  EXPECT_NEAR(-0.0141565285883 / sqrt(1.1930904), hessian.xz(2, 0), expectedPrecision);
  EXPECT_NEAR(-0.000302401003907 / sqrt(1.1930904), hessian.xz(3, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(1, 1), expectedPrecision);
  EXPECT_NEAR(0.000809443998852 / sqrt(1.1826193), hessian.xz(2, 1), expectedPrecision);
  EXPECT_NEAR(-0.105278309271 / sqrt(1.1826193), hessian.xz(3, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(1, 2), expectedPrecision);
  EXPECT_NEAR(-0.000285685983732 / sqrt(1.1826193), hessian.xz(2, 2), expectedPrecision);
  EXPECT_NEAR(-0.0300688478423 / sqrt(1.1826193), hessian.xz(3, 2), expectedPrecision);
  EXPECT_NEAR(1.23228518779 / sqrt(1.1826193), hessian.xz(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(1, 3), expectedPrecision);
  EXPECT_NEAR(-0.000142844126695 / sqrt(1.1826193), hessian.xz(2, 3), expectedPrecision);
  EXPECT_NEAR(0.0301435314879 / sqrt(1.1826193), hessian.xz(3, 3), expectedPrecision);
  EXPECT_NEAR(0.00628387012988 / sqrt(1.2623580), hessian.xz(0, 4), expectedPrecision);
  EXPECT_NEAR(-0.00569272697329 / sqrt(1.2623580), hessian.xz(1, 4), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(2, 4), expectedPrecision);
  EXPECT_NEAR(-0.000269390023198 / sqrt(1.2623580), hessian.xz(3, 4), expectedPrecision);
  EXPECT_NEAR(0.000661243126546 / sqrt(1.2623580), hessian.xz(0, 5), expectedPrecision);
  EXPECT_NEAR(0.0184515379115 / sqrt(1.2623580), hessian.xz(1, 5), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(2, 5), expectedPrecision);
  EXPECT_NEAR(0.000249962782101 / sqrt(1.2623580), hessian.xz(3, 5), expectedPrecision);
  EXPECT_NEAR(0.000154290062861 / sqrt(1.2623580), hessian.xz(0, 6), expectedPrecision);
  EXPECT_NEAR(0.0000258969272459 / sqrt(1.2623580), hessian.xz(1, 6), expectedPrecision);
  EXPECT_NEAR(-0.0444607527913 / sqrt(1.2623580), hessian.xz(2, 6), expectedPrecision);
  EXPECT_NEAR(0.000187472112098 / sqrt(1.2623580), hessian.xz(3, 6), expectedPrecision);
  EXPECT_NEAR(-0.0411830327973 / sqrt(1.2623580), hessian.xz(0, 7), expectedPrecision);
  EXPECT_NEAR(-0.0177376150546 / sqrt(1.2623580), hessian.xz(1, 7), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(2, 7), expectedPrecision);
  EXPECT_NEAR(-0.0000769685737339 / sqrt(1.2623580), hessian.xz(3, 7), expectedPrecision);
  EXPECT_NEAR(-0.000661243126546 / sqrt(1.2623580), hessian.xz(0, 8), expectedPrecision);
  EXPECT_NEAR(-0.0000361352473199 / sqrt(1.2623580), hessian.xz(1, 8), expectedPrecision);
  EXPECT_NEAR(-0.314599681951 / sqrt(1.2623580), hessian.xz(2, 8), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.xz(3, 8), expectedPrecision);
  EXPECT_NEAR(0.00628387012988 / sqrt(1.2623580), hessian.xz(0, 9), expectedPrecision);
  EXPECT_NEAR(0.00282132825293 / sqrt(1.2623580), hessian.xz(1, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.xz(2, 9), expectedPrecision);
  EXPECT_NEAR(-0.0000384843019961 / sqrt(1.2623580), hessian.xz(3, 9), expectedPrecision);
  EXPECT_NEAR(-1.90486012083 / sqrt(1.1930904), hessian.yy(0, 0), expectedPrecision);
  EXPECT_NEAR(-0.589265770897 / sqrt(1.1930904), hessian.yy(1, 0), expectedPrecision);
  EXPECT_NEAR(0.00354207798465 / sqrt(1.1930904), hessian.yy(2, 0), expectedPrecision);
  EXPECT_NEAR(-0.0000251994842742 / sqrt(1.1930904), hessian.yy(3, 0), expectedPrecision);
  EXPECT_NEAR(1.23228518779 / sqrt(1.1826193), hessian.yy(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.000142844126695 / sqrt(1.1826193), hessian.yy(2, 1), expectedPrecision);
  EXPECT_NEAR(-0.0149597402755 / sqrt(1.1826193), hessian.yy(3, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(1, 2), expectedPrecision);
  EXPECT_NEAR(-0.0000476131957936 / sqrt(1.1826193), hessian.yy(2, 2), expectedPrecision);
  EXPECT_NEAR(-0.0187370171672 / sqrt(1.1826193), hessian.yy(3, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(1, 3), expectedPrecision);
  EXPECT_NEAR(0.0000476147088982 / sqrt(1.1826193), hessian.yy(2, 3), expectedPrecision);
  EXPECT_NEAR(0.00747987013777 / sqrt(1.1826193), hessian.yy(3, 3), expectedPrecision);
  EXPECT_NEAR(0.0347464551086 / sqrt(1.2623580), hessian.yy(0, 4), expectedPrecision);
  EXPECT_NEAR(0.0298158834476 / sqrt(1.2623580), hessian.yy(1, 4), expectedPrecision);
  EXPECT_NEAR(0.289976751550 / sqrt(1.2623580), hessian.yy(2, 4), expectedPrecision);
  EXPECT_NEAR(0.000115452875730 / sqrt(1.2623580), hessian.yy(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.0947141661588 / sqrt(1.2623580), hessian.yy(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0400475913656 / sqrt(1.2623580), hessian.yy(1, 5), expectedPrecision);
  EXPECT_NEAR(0.124182336147 / sqrt(1.2623580), hessian.yy(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000333284093800 / sqrt(1.2623580), hessian.yy(3, 5), expectedPrecision);
  EXPECT_NEAR(-0.0601826256310 / sqrt(1.2623580), hessian.yy(0, 6), expectedPrecision);
  EXPECT_NEAR(-0.0258213125019 / sqrt(1.2623580), hessian.yy(1, 6), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.yy(3, 6), expectedPrecision);
  EXPECT_NEAR(0.0689968012890 / sqrt(1.2623580), hessian.yy(0, 7), expectedPrecision);
  EXPECT_NEAR(0.0140425974665 / sqrt(1.2623580), hessian.yy(1, 7), expectedPrecision);
  EXPECT_NEAR(-0.246517482483 / sqrt(1.2623580), hessian.yy(2, 7), expectedPrecision);
  EXPECT_NEAR(-4.81051962604e-6 / sqrt(1.2623580), hessian.yy(3, 7), expectedPrecision);
  EXPECT_NEAR(0.0947141661588 / sqrt(1.2623580), hessian.yy(0, 8), expectedPrecision);
  EXPECT_NEAR(0.0200237956828 / sqrt(1.2623580), hessian.yy(1, 8), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000166642046900 / sqrt(1.2623580), hessian.yy(3, 8), expectedPrecision);
  EXPECT_NEAR(0.0347464551086 / sqrt(1.2623580), hessian.yy(0, 9), expectedPrecision);
  EXPECT_NEAR(0.00745397086190 / sqrt(1.2623580), hessian.yy(1, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.yy(2, 9), expectedPrecision);
  EXPECT_NEAR(0.0000288632189325 / sqrt(1.2623580), hessian.yy(3, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(0, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(1, 0), expectedPrecision);
  EXPECT_NEAR(0.00707826429416 / sqrt(1.1930904), hessian.yz(2, 0), expectedPrecision);
  EXPECT_NEAR(-0.0000504001673179 / sqrt(1.1930904), hessian.yz(3, 0), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.000285685983732 / sqrt(1.1826193), hessian.yz(2, 1), expectedPrecision);
  EXPECT_NEAR(-0.0300688478423 / sqrt(1.1826193), hessian.yz(3, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(1, 2), expectedPrecision);
  EXPECT_NEAR(0.0000476147088982 / sqrt(1.1826193), hessian.yz(2, 2), expectedPrecision);
  EXPECT_NEAR(0.00747987013777 / sqrt(1.1826193), hessian.yz(3, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(1, 3), expectedPrecision);
  EXPECT_NEAR(0.0000476147088982 / sqrt(1.1826193), hessian.yz(2, 3), expectedPrecision);
  EXPECT_NEAR(0.00753588287198 / sqrt(1.1826193), hessian.yz(3, 3), expectedPrecision);
  EXPECT_NEAR(0.0164732131189 / sqrt(1.2623580), hessian.yz(0, 4), expectedPrecision);
  EXPECT_NEAR(0.0141900920437 / sqrt(1.2623580), hessian.yz(1, 4), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(2, 4), expectedPrecision);
  EXPECT_NEAR(0.000153937147468 / sqrt(1.2623580), hessian.yz(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.0601826256310 / sqrt(1.2623580), hessian.yz(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0258213125019 / sqrt(1.2623580), hessian.yz(1, 5), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.yz(3, 5), expectedPrecision);
  EXPECT_NEAR(-0.000661243126546 / sqrt(1.2623580), hessian.yz(0, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000361352473199 / sqrt(1.2623580), hessian.yz(1, 6), expectedPrecision);
  EXPECT_NEAR(-0.314599681951 / sqrt(1.2623580), hessian.yz(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.yz(3, 6), expectedPrecision);
  EXPECT_NEAR(0.0707746935497 / sqrt(1.2623580), hessian.yz(0, 7), expectedPrecision);
  EXPECT_NEAR(0.0150978354912 / sqrt(1.2623580), hessian.yz(1, 7), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(2, 7), expectedPrecision);
  EXPECT_NEAR(0.0000192421509981 / sqrt(1.2623580), hessian.yz(3, 7), expectedPrecision);
  EXPECT_NEAR(0.00154290062861 / sqrt(1.2623580), hessian.yz(0, 8), expectedPrecision);
  EXPECT_NEAR(0.0000421577885398 / sqrt(1.2623580), hessian.yz(1, 8), expectedPrecision);
  EXPECT_NEAR(-0.306627154417 / sqrt(1.2623580), hessian.yz(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000374944271145 / sqrt(1.2623580), hessian.yz(3, 8), expectedPrecision);
  EXPECT_NEAR(-0.0157096753247 / sqrt(1.2623580), hessian.yz(0, 9), expectedPrecision);
  EXPECT_NEAR(-0.00352666031616 / sqrt(1.2623580), hessian.yz(1, 9), expectedPrecision);
  EXPECT_NEAR(0, hessian.yz(2, 9), expectedPrecision);
  EXPECT_NEAR(0.0000192421509981 / sqrt(1.2623580), hessian.yz(3, 9), expectedPrecision);
  EXPECT_NEAR(-1.90486012083 / sqrt(1.1930904), hessian.zz(0, 0), expectedPrecision);
  EXPECT_NEAR(-0.589265770897 / sqrt(1.1930904), hessian.zz(1, 0), expectedPrecision);
  EXPECT_NEAR(0.00354207798465 / sqrt(1.1930904), hessian.zz(2, 0), expectedPrecision);
  EXPECT_NEAR(0.0000504007667026 / sqrt(1.1930904), hessian.zz(3, 0), expectedPrecision);
  EXPECT_NEAR(1.23228518779 / sqrt(1.1826193), hessian.zz(0, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(1, 1), expectedPrecision);
  EXPECT_NEAR(-0.000142844126695 / sqrt(1.1826193), hessian.zz(2, 1), expectedPrecision);
  EXPECT_NEAR(0.0301435314879 / sqrt(1.1826193), hessian.zz(3, 1), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(0, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(1, 2), expectedPrecision);
  EXPECT_NEAR(0.0000476147088982 / sqrt(1.1826193), hessian.zz(2, 2), expectedPrecision);
  EXPECT_NEAR(0.00753588287198 / sqrt(1.1826193), hessian.zz(3, 2), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(0, 3), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(1, 3), expectedPrecision);
  EXPECT_NEAR(-0.0000476131957936 / sqrt(1.1826193), hessian.zz(2, 3), expectedPrecision);
  EXPECT_NEAR(0.0149223984527 / sqrt(1.1826193), hessian.zz(3, 3), expectedPrecision);
  EXPECT_NEAR(0.000152707558845 / sqrt(1.2623580), hessian.zz(0, 4), expectedPrecision);
  EXPECT_NEAR(0.0000166901558139 / sqrt(1.2623580), hessian.zz(1, 4), expectedPrecision);
  EXPECT_NEAR(-0.121089474041 / sqrt(1.2623580), hessian.zz(2, 4), expectedPrecision);
  EXPECT_NEAR(0.000115452875730 / sqrt(1.2623580), hessian.zz(3, 4), expectedPrecision);
  EXPECT_NEAR(-0.000661243126546 / sqrt(1.2623580), hessian.zz(0, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000361352473199 / sqrt(1.2623580), hessian.zz(1, 5), expectedPrecision);
  EXPECT_NEAR(-0.314599681951 / sqrt(1.2623580), hessian.zz(2, 5), expectedPrecision);
  EXPECT_NEAR(-0.0000999851233221 / sqrt(1.2623580), hessian.zz(3, 5), expectedPrecision);
  EXPECT_NEAR(0.0220324619169 / sqrt(1.2623580), hessian.zz(0, 6), expectedPrecision);
  EXPECT_NEAR(0.00978782185672 / sqrt(1.2623580), hessian.zz(1, 6), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(2, 6), expectedPrecision);
  EXPECT_NEAR(-0.0000333284093800 / sqrt(1.2623580), hessian.zz(3, 6), expectedPrecision);
  EXPECT_NEAR(0.000954422242778 / sqrt(1.2623580), hessian.zz(0, 7), expectedPrecision);
  EXPECT_NEAR(0.0000260783684592 / sqrt(1.2623580), hessian.zz(1, 7), expectedPrecision);
  EXPECT_NEAR(-0.272451316592 / sqrt(1.2623580), hessian.zz(2, 7), expectedPrecision);
  EXPECT_NEAR(0.0000288632189325 / sqrt(1.2623580), hessian.zz(3, 7), expectedPrecision);
  EXPECT_NEAR(-0.0550811547921 / sqrt(1.2623580), hessian.zz(0, 8), expectedPrecision);
  EXPECT_NEAR(-0.0122347773209 / sqrt(1.2623580), hessian.zz(1, 8), expectedPrecision);
  EXPECT_NEAR(0, hessian.zz(2, 8), expectedPrecision);
  EXPECT_NEAR(0.0000166642046900 / sqrt(1.2623580), hessian.zz(3, 8), expectedPrecision);
  EXPECT_NEAR(-0.0128477041176 / sqrt(1.2623580), hessian.zz(0, 9), expectedPrecision);
  EXPECT_NEAR(-0.00282897790767 / sqrt(1.2623580), hessian.zz(1, 9), expectedPrecision);
  EXPECT_NEAR(0.190840092897 / sqrt(1.2623580), hessian.zz(2, 9), expectedPrecision);
  EXPECT_NEAR(-4.81051962604e-6 / sqrt(1.2623580), hessian.zz(3, 9), expectedPrecision);
}

/* =========================
 *   Spherical Basis Tests
 * =========================
 *  ( Minimal test of normalization/completeness. June 06, 2017 - JU)
 */

class SphBFOnGridTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    std::string baslibpath;
    if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
      baslibpath = (std::string)env_p + "basis/";
    }
    else {
      throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
    }

    _geom = std::make_shared<Geometry>(std::vector<std::shared_ptr<Atom>>(1, std::make_shared<Atom>("SC", 0.0, 0.0, 0.0)));
    auto basis = AtomCenteredBasisControllerFactory::produce(_geom, baslibpath, true, true, 37, "CC-PV5Z");

    auto grid = AtomCenteredGridControllerFactory::produce(_geom, Options::GRID_TYPES::BECKE, 3,
                                                           Options::RADIAL_GRID_TYPES::AHLRICHS,
                                                           Options::SPHERICAL_GRID_TYPES::LEBEDEV, 7, 0.0, true); /*1e-14*/

    _weights = grid->getWeights();

    auto bfOnGridCalc = BasisFunctionOnGridControllerFactory::produce(33310, 0.0, 1, basis, grid); /* 1e-09 */
    _data = bfOnGridCalc->getBlockOnGridData(0);
    std::remove("WARNING");
  }
  static std::shared_ptr<Geometry> _geom;
  static Eigen::VectorXd _weights;
  static std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> _data;
  double _expectedPrecision = 1E-7;
};
std::shared_ptr<Geometry> SphBFOnGridTest::_geom = nullptr;
Eigen::VectorXd SphBFOnGridTest::_weights = Eigen::VectorXd();
std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> SphBFOnGridTest::_data = nullptr;
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, S_Functions) {
  for (unsigned int n = 0; n < 9; n++) { // for n>7 the grid is not large enough
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, S_Functions_Grad) {
  for (unsigned int n = 0; n < 9; n++) { // gradients are 0 for s functions so they can run till n=9
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, P_Functions) {
  for (unsigned int n = 9; n < 33; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, P_Functions_Grad) {
  for (unsigned int n = 9; n < 30; n++) { // for n>31 the grid is not large enough
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision * 100.0);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision * 100.0);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision * 100.0);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, D_Functions) {
  for (unsigned int n = 33; n < 63; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision * 20.0);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, D_Functions_Grad) {
  for (unsigned int n = 33; n < 63; n++) {
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision * 10.0);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision * 10.0);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision * 10.0);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, F_Functions) {
  for (unsigned int n = 63; n < 91; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, F_Functions_Grad) {
  for (unsigned int n = 63; n < 91; n++) {
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, G_Functions) {
  for (unsigned int n = 91; n < 118; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, G_Functions_Grad) {
  for (unsigned int n = 91; n < 118; n++) {
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, H_Functions) {
  for (unsigned int n = 118; n < 140; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, H_Functions_Grad) {
  for (unsigned int n = 118; n < 140; n++) {
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, I_Functions) {
  for (unsigned int n = 140; n < 153; n++) {
    double integral_f_r =
        (_data->functionValues.col(n).array() * _data->functionValues.col(n).array() * _weights.array()).sum();
    EXPECT_NEAR(1.0, integral_f_r, _expectedPrecision);
  }
}
/**
 * @test
 * @brief Minimal test of normalization/completeness.
 */
TEST_F(SphBFOnGridTest, I_Functions_Grad) {
  for (unsigned int n = 140; n < 153; n++) {
    EXPECT_NEAR(0.0, _data->derivativeValues->x.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->y.col(n).dot(_weights), _expectedPrecision);
    EXPECT_NEAR(0.0, _data->derivativeValues->z.col(n).dot(_weights), _expectedPrecision);
  }
}
/**
 * @test
 * @brief See that 2nd derirvs are possible and do not crash.
 */
TEST_F(SphBFOnGridTest, Hess_Calc) {
  std::string baslibpath;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    baslibpath = (std::string)env_p + "basis/";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }

  auto basis = AtomCenteredBasisControllerFactory::produce(_geom, baslibpath, true, true, 37, "CC-PV5Z");

  auto grid = AtomCenteredGridControllerFactory::produce(_geom, Options::GRID_TYPES::BECKE, 3,
                                                         Options::RADIAL_GRID_TYPES::AHLRICHS,
                                                         Options::SPHERICAL_GRID_TYPES::LEBEDEV, 4, 0.0, true); /*1e-14*/

  auto bfOnGridCalc = BasisFunctionOnGridControllerFactory::produce(128, 0.0, 2, basis, grid); /* 1e-09 */
  bfOnGridCalc->getBlockOnGridData(0);
  bfOnGridCalc->getBlockOnGridData(1);
}
} /* namespace Serenity */
