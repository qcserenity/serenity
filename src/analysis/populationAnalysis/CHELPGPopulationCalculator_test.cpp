/**
 * @file   CHELPGPopulationCalculator_test.cpp
 *
 * @date   Oct 7, 2024
 * @author Thorben Wiegmann
 *
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
#include "analysis/populationAnalysis/CHELPGPopulationCalculator.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory>

namespace Serenity {

/**
 * @class CHELPGPopulationCalculatorTest
 * @brief Sets everything up for the tests of CHELPGPopulationCalculator.h/.cpp.\n
          All of these results were obtained with ORCA 6.0.0 on October 07, 2024.
 */
class CHELPGPopulationCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(CHELPGPopulationCalculatorTest, waterHFRestricted) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  CHELPGPopulationCalculator<Options::SCF_MODES::RESTRICTED> calculator(system);
  auto populations = calculator.getAtomPopulations();
  /*
   * Check results
   */
  EXPECT_NEAR(populations[0], 0.401188, 1e-2);
  EXPECT_NEAR(populations[1], 0.401425, 1e-2);
  EXPECT_NEAR(populations[2], -0.802614, 1e-2);
}

TEST_F(CHELPGPopulationCalculatorTest, waterHFUnrestricted) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  CHELPGPopulationCalculator<Options::SCF_MODES::UNRESTRICTED> calculator(system);
  auto populations = calculator.getAtomPopulations();
  /*
   * Check results
   */
  EXPECT_NEAR(populations[0], 0.401215, 1e-2);
  EXPECT_NEAR(populations[1], 0.401452, 1e-2);
  EXPECT_NEAR(populations[2], -0.802667, 1e-2);
}

TEST_F(CHELPGPopulationCalculatorTest, waterDFTRestricted) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  CHELPGPopulationCalculator<Options::SCF_MODES::RESTRICTED> calculator(system);
  auto populations = calculator.getAtomPopulations();
  /*
   * Check results
   */
  EXPECT_NEAR(populations[0], 0.371879, 1e-2);
  EXPECT_NEAR(populations[1], 0.372103, 1e-2);
  EXPECT_NEAR(populations[2], -0.743979, 1e-2);
}

TEST_F(CHELPGPopulationCalculatorTest, waterDFTUnrestricted) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  CHELPGPopulationCalculator<Options::SCF_MODES::UNRESTRICTED> calculator(system);
  auto populations = calculator.getAtomPopulations();
  /*
   * Check results
   */
  EXPECT_NEAR(populations[0], 0.371876, 1e-2);
  EXPECT_NEAR(populations[1], 0.372103, 1e-2);
  EXPECT_NEAR(populations[2], -0.743979, 1e-2);
}

} // namespace Serenity
