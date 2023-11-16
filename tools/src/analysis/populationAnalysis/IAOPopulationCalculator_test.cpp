/**
 * @file   IAOPopulationCalculator_test.cpp
 *
 * @date   Oct 17, 2018
 * @author Moritz Bensberg
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
#include "analysis/populationAnalysis/IAOPopulationCalculator.h" //To be tested.
#include "testsupply/SystemController__TEST_SUPPLY.h"            //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <memory>        //smrt_ptr.

namespace Serenity {

class IAOPopulationCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Checks correct functionality of the atom population calculator for a water molecule.
 */
TEST_F(IAOPopulationCalculatorTest, unrestricted_atomPopulations) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);

  auto populations = IAOPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>::calculateIAOPopulations(system);
  /*
   * Check results
   */
  EXPECT_NEAR(populations.alpha[0], 0.65211998967722806 / 2.0, 1e-5);
  EXPECT_NEAR(populations.alpha[1], 0.65212009368996549 / 2.0, 1e-5);
  EXPECT_NEAR(populations.alpha[2], 8.695759916632829 / 2.0, 1e-5);
}
/**
 * @test
 * @brief Checks correct functionality of the atom population calculator for a water molecule.
 */
TEST_F(IAOPopulationCalculatorTest, restricted_atomPopulations) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);

  auto populations = IAOPopulationCalculator<Options::SCF_MODES::RESTRICTED>::calculateIAOPopulations(system);
  /*
   * Check results
   */
  EXPECT_NEAR(populations[0], 0.65211998967722806, 1e-5);
  EXPECT_NEAR(populations[1], 0.65212009368996549, 1e-5);
  EXPECT_NEAR(populations[2], 8.695759916632829, 1e-5);
}

/**
 * @test
 * @brief Checks correct functionality of the atom population calculator for a water molecule.
 */
TEST_F(IAOPopulationCalculatorTest, restricted_1SPopulations) {
  /*
   * Set up test environment
   */
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);

  auto populations = IAOPopulationCalculator<Options::SCF_MODES::RESTRICTED>::calculate1SOrbitalPopulations(system);
  /*
   * Check results
   */
  EXPECT_NEAR(populations.col(0).array().sum(), 1.0, 1e-3); // 1s on oxygen.
  EXPECT_NEAR(populations.col(4).array().sum(), 0.0, 1e-3); // Lone pair on oxygen.
}
} /* namespace Serenity */
