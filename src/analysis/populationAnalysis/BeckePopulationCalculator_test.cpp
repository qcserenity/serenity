/**
 * @file   BeckePopulationCalculator_test.cpp
 *
 * @date   Sep 15, 2020
 * @author Patrick Eschenbach
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
#include "analysis/populationAnalysis/BeckePopulationCalculator.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class BeckePopulationCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Checks atom density populations for restricted CO minbas
 */
TEST_F(BeckePopulationCalculatorTest, populationsRestricted) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto densPtr = std::make_shared<DensityMatrix<RESTRICTED>>(
      systemController->getElectronicStructure<RESTRICTED>()->getDensityMatrix());
  auto beckeCalculator = BeckePopulationCalculator<RESTRICTED>(systemController, densPtr);
  auto atomPops = beckeCalculator.getAtomPopulations();
  EXPECT_NEAR(atomPops[0], 5.73331, 1e-4);
  EXPECT_NEAR(atomPops[1], 8.26668, 1e-4);
}

/**
 * @test
 * @brief Checks atom spin populations for H2_DEF2_TZVP_PW91_UNRES_1
 */
TEST_F(BeckePopulationCalculatorTest, spinPopulationsUnrestricted) {
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);
  auto densPtr = std::make_shared<DensityMatrix<UNRESTRICTED>>(
      systemController->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix());
  auto beckeCalculator = BeckePopulationCalculator<UNRESTRICTED>(systemController, densPtr);
  auto spinPops = beckeCalculator.getSpinPopulations();
  EXPECT_NEAR(spinPops[0], 0.5, 1e-6);
  EXPECT_NEAR(spinPops[1], 0.5, 1e-6);
}

} /* namespace Serenity */
