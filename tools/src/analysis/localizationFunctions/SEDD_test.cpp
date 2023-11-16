/**
 * @file SEDD_test.cpp
 *
 * @date Jul 12, 2017
 * @author Jan Unsleber
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
#include "analysis/localizationFunctions/SEDD.h"
#include "data/grid/GridData.h"
#include "grid/GridController.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class SEDDTest : public ::testing::Test {
 protected:
  SEDDTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS)) {
    grid = systemController->getGridController(Options::GRID_PURPOSES::SMALL);
  }
  virtual ~SEDDTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
  // system
  std::shared_ptr<SystemController> systemController;
  // grid
  std::shared_ptr<GridController> grid;
};

/**
 * @test
 * @brief Tests the calculation of the SEDD of CO in a minimal basis
 */
TEST_F(SEDDTest, SEDD_restricted) {
  SEDD<Options::SCF_MODES::RESTRICTED> sedd;
  auto lambda = sedd.getSEDDLambda(systemController);
  auto result = lambda(grid);
  EXPECT_NEAR(result.sum(), 55642.7514261, 5e-6);
}
/**
 * @test
 * @brief Tests the calculation of the SEDD of CO in a minimal basis
 */
TEST_F(SEDDTest, SEDD_unrestricted) {
  SEDD<Options::SCF_MODES::UNRESTRICTED> sedd;
  auto lambda = sedd.getSEDDLambda(systemController);
  auto result = lambda(grid);
  EXPECT_NEAR(result.sum(), 55642.7514261, 5e-6);
}
/**
 * @test
 * @brief Tests the calculation of the DORI of CO in a minimal basis
 */
TEST_F(SEDDTest, DORI_restricted) {
  SEDD<Options::SCF_MODES::RESTRICTED> dori;
  auto lambda = dori.getDORILambda(systemController);
  auto result = lambda(grid);
  EXPECT_TRUE(result.minCoeff() >= 0.0);
  EXPECT_NEAR(result.sum(), 1428.4446333128767, 5e-6);
}
/**
 * @test
 * @brief Tests the calculation of the DORI of CO in a minimal basis
 */
TEST_F(SEDDTest, DORI_unrestricted) {
  SEDD<Options::SCF_MODES::UNRESTRICTED> dori;
  auto lambda = dori.getDORILambda(systemController);
  auto result = lambda(grid);
  EXPECT_TRUE(result.minCoeff() >= 0.0);
  EXPECT_NEAR(result.sum(), 1428.4446333128767, 5e-6);
}
/**
 * @test
 * @brief Tests the calculation of the signed density of CO in a minimal basis
 */
TEST_F(SEDDTest, SignedDensity_restricted) {
  SEDD<Options::SCF_MODES::RESTRICTED> dori;
  auto lambda = dori.getSignedDensityLambda(systemController);
  auto result = lambda(grid);
  EXPECT_NEAR(result.sum(), -106849.88204633593, 5e-6);
}
/**
 * @test
 * @brief Tests the calculation of the the signed density of CO in a minimal basis
 */
TEST_F(SEDDTest, SignedDensity_unrestricted) {
  SEDD<Options::SCF_MODES::UNRESTRICTED> dori;
  auto lambda = dori.getSignedDensityLambda(systemController);
  auto result = lambda(grid);
  EXPECT_NEAR(result.sum(), -106849.88204633593, 5e-6);
}

} // namespace Serenity
