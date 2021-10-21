/**
 * @file CC_test.cpp
 *
 * @date Apr 14, 2016
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
#include "data/ElectronicStructure.h"
#include "postHF/CC/CCSD_T.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class CoupledCluster : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Tests the R-CCSD(T) result for distorted water in a minimal basis.
 */
TEST_F(CoupledCluster, RCCSD_T_Water_MinBas) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  CCSD_T ccsd_t(systemController, 1e-7, 100);

  auto result = ccsd_t.calculateElectronicEnergyCorrections();
  auto triples = ccsd_t.calculateTripplesCorrection();
  // check MP2 energy
  EXPECT_NEAR(result.first, -0.027671561078791318, 1e-6);
  // check R-CCSD energy
  EXPECT_NEAR(result.second, -0.038242629623747561, 1e-6);
  // check R-CCSD(T) tripples
  EXPECT_NEAR(triples, -0.00016316279608078255, 1e-6);
}

/**
 * @test
 * @brief Tests the R-CCSD result for H2 in a minimal basis.
 */
TEST_F(CoupledCluster, RCCSD_H2_MinimalBasis) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  CCSD ccsd(systemController, 1e-7, 100);

  auto result = ccsd.calculateElectronicEnergyCorrections();

  // check MP2 energy
  EXPECT_NEAR(result.first, -0.01407063069504663, 1e-6);
  // check R-CCSD energy
  EXPECT_NEAR(result.second, -0.022809684727466688, 1e-6);
}

} /* namespace Serenity */
