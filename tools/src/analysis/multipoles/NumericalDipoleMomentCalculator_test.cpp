/**
 * @file NumericalDipoleMomentCalculator_test.cpp
 *
 * @date Oct 30, 2015
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
#include "analysis/multipoles/NumericalDipoleMomentCalculator.h"
#include "settings/Options.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class NumericalDipoleMomentCalculatorTest : public ::testing::Test {
 protected:
  NumericalDipoleMomentCalculator calc = NumericalDipoleMomentCalculator();
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Numerical calculation of the dipole moment of H2 based on the density on a grid
 *
 * unrestricted, minimal basis
 */
TEST_F(NumericalDipoleMomentCalculatorTest, H2_Unrestricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto results = calc.calculateDipoleMoment<Options::SCF_MODES::UNRESTRICTED>(system);
  /*
   * Check results
   */
  EXPECT_EQ((unsigned int)3, results.size());
  EXPECT_NEAR(-0.000000000, results[0], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[1], 1.0e-6);
  EXPECT_NEAR(+0.000000000, results[2], 1.0e-6);
};

/**
 * @test
 * @brief Numerical calculation of the dipole moment of H2 based on the density on a grid
 *
 * restricted, minimal basis
 */
TEST_F(NumericalDipoleMomentCalculatorTest, H2_Restricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto results = calc.calculateDipoleMoment<Options::SCF_MODES::RESTRICTED>(system);
  /*
   * Check results
   */
  EXPECT_EQ((unsigned int)3, results.size());
  EXPECT_NEAR(-0.000000000, results[0], 1.0e-6);
  EXPECT_NEAR(-0.000000000, results[1], 1.0e-6);
  EXPECT_NEAR(+0.000000000, results[2], 1.0e-6);
};

/**
 * @test
 * @brief Numerical calculation of the dipole moment of CO based on the density on a grid
 *
 * minimal basis
 */
TEST_F(NumericalDipoleMomentCalculatorTest, CO) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto results = calc.calculateDipoleMoment<Options::SCF_MODES::RESTRICTED>(system);
  /*
   * Check results
   */
  EXPECT_EQ((unsigned int)3, results.size());
  EXPECT_NEAR(-0.000000000, results[0], 1.0e-5);
  EXPECT_NEAR(-0.000000000, results[1], 1.0e-5);
  EXPECT_NEAR(+0.040437981, results[2], 1.0e-5);
}

} /* namespace Serenity */
