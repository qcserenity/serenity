/**
 * @file LocalCorrelationTask_test.cpp
 *
 * @date Apr 4, 2021
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
#include "tasks/LocalCorrelationTask.h"               //To be tested.
#include "data/ElectronicStructure.h"                 //Access to energies.
#include "energies/EnergyContributions.h"             //Access to energies.
#include "system/SystemController.h"                  //Access to electronic structure.
#include "tasks/ScfTask.h"                            //Run SCF and ignore previous electronic structure.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LocalCorrelationTaskTest
 * @brief Sets everything up for the tests of LocalCorrelationTask.h/.cpp .
 */
class LocalCorrelationTaskTest : public ::testing::Test {
 protected:
  LocalCorrelationTaskTest() {
  }

  virtual ~LocalCorrelationTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Test DLPNO-CCSD(T0) with tight settings.
 *        Note, this should be identical to the test CoupledClusterTaskTest.WaterDLPNOCCSDT0_TIGHT.
 */
TEST_F(LocalCorrelationTaskTest, water_ccsdt0_tight) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalCorrelationTask task(act);
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(-0.215243, ccsdCorrection, 1e-6);
  // Canonical: -0.00329571. However, this is only semi-canonical treatment!
  EXPECT_NEAR(-0.003178, triplesCorrection, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Test DLPNO-CCSD(T0) with tight settings.
 *        Note, this should be identical to the test CoupledClusterTaskTest.WaterDLPNOCCSDT0_LOOSE.
 */
TEST_F(LocalCorrelationTaskTest, water_ccsdt0_loose) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalCorrelationTask task(act);
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(-0.215378, ccsdCorrection, 1e-6);
  // Canonical: -0.00329571. However, this is only semi-canonical treatment!
  EXPECT_NEAR(-0.003153, triplesCorrection, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Test DLPNO-MP2 with tight settings.
 *        Note, this should be identical to the test MP2TaskTest.LocalMP2Restricted_TIGHT.
 */
TEST_F(LocalCorrelationTaskTest, water_mp2_tight) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  LocalCorrelationTask task(act);
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  task.settings.lcSettings.reuseFockMatrix = false;
  task.run();
  auto energyComponentController =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.205918, 1E-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Test DLPNO-CCSD(T0) with tight settings. Check if the triples correction gives zero.
 */
TEST_F(LocalCorrelationTaskTest, h2_ccsdt0_tight) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalCorrelationTask task(act);
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  // Canonical CCSD reference.
  EXPECT_NEAR(-0.035726, ccsdCorrection, 1e-4);
  EXPECT_NEAR(0.0, triplesCorrection, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/