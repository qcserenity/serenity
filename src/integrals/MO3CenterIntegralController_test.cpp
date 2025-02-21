/**
 * @file MO3CenterIntegralController_test.cpp
 *
 * @date May 15, 2019
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
#include "integrals/MO3CenterIntegralController.h"              //To be tested.
#include "basis/BasisController.h"                              //Number of aux. basis functions.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Easy construction of the MO3CenterIntegralController.
#include "system/SystemController.h"                            //Test systems.
#include "testsupply/SystemController__TEST_SUPPLY.h"           //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Test framework.
namespace Serenity {
/**
 * @class MO3CenterIntegralControllerTest MO3CenterIntegralController_test.cpp
 * @brief Sets up the test of the MO3CenterIntegralController.
 *
 * Note: The class MO3CenterIntegralController is indirectly tested by the tests in LocalMP2Test, DLPNO_CCSDTest etc.
 */
class MO3CenterIntegralControllerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests the integrals against precalculated values with a second RI implementation.
 */
TEST_F(MO3CenterIntegralControllerTest, H2_RIvsRI) {
  // Build reference
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  auto aoBasis = system->getBasisController();

  LocalCorrelationSettings lCSettings;
  LocalCorrelationController lCController(system, lCSettings);

  auto auxBasis = system->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  auto m = auxBasis->getReducedNBasisFunctions();
  auto mo3CenterIntController = lCController.getMO3CenterIntegralController();
  auto exchangeInts =
      mo3CenterIntController->getMO3CenterInts(MO3CENTER_INTS::ia_K, Eigen::VectorXi::Constant(m, 1).sparseView());
  EXPECT_NEAR(exchangeInts[0](0, 0), -0.379003, 1e-6);
  EXPECT_NEAR(exchangeInts[0](1, 0), -0.0626244, 1e-6);
  EXPECT_NEAR(exchangeInts[0](2, 0), 0.181815, 1e-6);

  EXPECT_NEAR(exchangeInts[1](0, 0), -0.737497, 1e-6);
  EXPECT_NEAR(exchangeInts[1](1, 0), -0.200962, 1e-6);
  EXPECT_NEAR(exchangeInts[1](2, 0), 0.431401, 1e-6);

  EXPECT_NEAR(exchangeInts[2](0, 0), -0.784675, 1e-6);
  EXPECT_NEAR(exchangeInts[2](1, 0), -0.314367, 1e-6);
  EXPECT_NEAR(exchangeInts[2](2, 0), 0.581753, 1e-6);

  EXPECT_NEAR(exchangeInts[3](9, 0), 0.0, 1e-6);
  EXPECT_NEAR(exchangeInts[3](10, 0), 0.0, 1e-6);
  EXPECT_NEAR(exchangeInts[3](11, 0), 0.551752, 1e-6);
};

/**
 * @test
 * @brief Tests the writing and reading from disk.
 */
TEST_F(MO3CenterIntegralControllerTest, ReadAndWrite) {
  // Build reference
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  auto aoBasis = system->getBasisController();

  LocalCorrelationSettings lCSettings;
  LocalCorrelationController lCController(system, lCSettings);

  auto auxBasis = system->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  auto m = auxBasis->getReducedNBasisFunctions();
  auto mo3CenterIntController = lCController.getMO3CenterIntegralController();
  auto exchangeInts =
      mo3CenterIntController->getMO3CenterInts(MO3CENTER_INTS::ia_K, Eigen::VectorXi::Constant(m, 1).sparseView());
  mo3CenterIntController->removeFromMemory(MO3CENTER_INTS::ia_K);
  auto exchangeIntsFromDisk =
      mo3CenterIntController->getMO3CenterInts(MO3CENTER_INTS::ia_K, Eigen::VectorXi::Constant(m, 1).sparseView());
  EXPECT_NEAR(exchangeInts[0](0, 0), exchangeIntsFromDisk[0](0, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[0](1, 0), exchangeIntsFromDisk[0](1, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[0](2, 0), exchangeIntsFromDisk[0](2, 0), 1e-12);

  EXPECT_NEAR(exchangeInts[1](0, 0), exchangeIntsFromDisk[1](0, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[1](1, 0), exchangeIntsFromDisk[1](1, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[1](2, 0), exchangeIntsFromDisk[1](2, 0), 1e-12);

  EXPECT_NEAR(exchangeInts[2](0, 0), exchangeIntsFromDisk[2](0, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[2](1, 0), exchangeIntsFromDisk[2](1, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[2](2, 0), exchangeIntsFromDisk[2](2, 0), 1e-12);

  EXPECT_NEAR(exchangeInts[3](9, 0), exchangeIntsFromDisk[3](9, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[3](10, 0), exchangeIntsFromDisk[3](10, 0), 1e-12);
  EXPECT_NEAR(exchangeInts[3](11, 0), exchangeIntsFromDisk[3](11, 0), 1e-12);
};

} /* namespace Serenity */
