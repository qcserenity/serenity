/**
 * @file MultipoleMomentCalculator_test.cpp
 *
 * @date Dec 8, 2015
 * @author David Schnieders
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
#include "analysis/multipoles/MultipoleMomentCalculator.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class MultipoleMomentCalculatorTest : public ::testing::Test {
 protected:
  MultipoleMomentCalculator calc = MultipoleMomentCalculator();

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Analytical calculation of the multipole moments of H2 based on the density on a grid
 *
 * unrestricted, minimal basis, uses old Serenity results as reference values
 */
TEST_F(MultipoleMomentCalculatorTest, H2_Unrestricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto results = calc.calculateMultipoleMoment<Options::SCF_MODES::UNRESTRICTED>(system, 2);
  /*
   * Check results
   */
  EXPECT_NEAR(-0.000000000, results[0][0], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[0][1], 1.0e-6);
  EXPECT_NEAR(+0.000000000, results[0][2], 1.0e-6);
  EXPECT_NEAR(-1.37665417, results[1][0], 1.0e-6);
  EXPECT_NEAR(-1.06516997, results[1][5], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[1][1], 1.0e-6);
};

/**
 * @test
 * @brief Analytical calculation of the multipole moments of H2 based on the density on a grid
 *
 * restricted, minimal basis, uses old Serenity results as reference values
 */
TEST_F(MultipoleMomentCalculatorTest, H2_Restricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto results = calc.calculateMultipoleMoment<Options::SCF_MODES::RESTRICTED>(system, 2);
  /*
   * Check results
   */
  EXPECT_NEAR(-0.000000000, results[0][0], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[0][1], 1.0e-6);
  EXPECT_NEAR(+0.000000000, results[0][2], 1.0e-6);
  EXPECT_NEAR(-1.37665417, results[1][0], 1.0e-6);
  EXPECT_NEAR(-1.06516997, results[1][5], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[1][1], 1.0e-6);
};

/**
 * @test
 * @brief Analytical calculation of the multipole moments of H2 based on the density on a grid
 *
 * minimal basis, uses old Serenity results as reference values
 */
TEST_F(MultipoleMomentCalculatorTest, CO) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto results = calc.calculateMultipoleMoment<Options::SCF_MODES::RESTRICTED>(system, 2);
  /*
   * Check results
   */
  EXPECT_NEAR(-0.000000000, results[0][0], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[0][1], 1.0e-6);
  EXPECT_NEAR(+0.040437273, results[0][2], 1.0e-6);
  EXPECT_NEAR(-6.44718475, results[1][0], 1.0e-6);
  EXPECT_NEAR(-8.59555897, results[1][5], 1.0e-6);
  EXPECT_NEAR(+0.000000000, results[1][1], 1.0e-6);
}

} /* namespace Serenity */
