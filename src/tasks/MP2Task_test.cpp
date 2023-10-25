/**
 * @file MP2Task_test.cpp
 *
 * @date Aug 7, 2019
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
#include "tasks/MP2Task.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class MP2TaskTest
 * @brief Sets everything up for the tests of MP2Task.h/.cpp .
 */
class MP2TaskTest : public ::testing::Test {
 protected:
  MP2TaskTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS)) {
    systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  }

  virtual ~MP2TaskTest() = default;

  /// system
  std::shared_ptr<SystemController> systemController;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test.
 */
TEST_F(MP2TaskTest, MP2Restricted) {
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(systemController, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::AO;
  mp2Task.run();
  auto energyComponentController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.01407063069504663, 1E-8);
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with RI.
 */
TEST_F(MP2TaskTest, MP2RestrictedRI) {
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(systemController, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::DF;
  mp2Task.run();
  auto energyComponentController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.014069347723441622, 1E-8);
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with LOCAL-MP2/LOOSE-PNO.
 */
TEST_F(MP2TaskTest, LocalMP2Restricted_LOOSE) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(system, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  mp2Task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  mp2Task.settings.lcSettings.reuseFockMatrix = false;
  mp2Task.run();
  auto energyComponentController =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  // Difference to full RI ~ 1e-4
  EXPECT_NEAR(MP2EnergyCorrection, -0.205835, 1E-6);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with LOCAL-MP2/NORMAL-PNO.
 */
TEST_F(MP2TaskTest, LocalMP2Restricted_NORMAL) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(system, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  mp2Task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::NORMAL;
  mp2Task.settings.lcSettings.reuseFockMatrix = false;
  mp2Task.run();
  auto energyComponentController =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  // Difference to full RI ~ 1e-4
  EXPECT_NEAR(MP2EnergyCorrection, -0.205866, 1E-6);
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with LOCAL-MP2/TIGHT-PNO.
 */
TEST_F(MP2TaskTest, LocalMP2Restricted_TIGHT) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(system, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  mp2Task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  mp2Task.settings.lcSettings.pnoThreshold = 1e-12;
  mp2Task.settings.lcSettings.orbitalToShellThreshold = 1e-5;
  mp2Task.settings.lcSettings.mullikenThreshold = 1e-5;
  mp2Task.settings.lcSettings.reuseFockMatrix = false;
  mp2Task.run();
  auto energyComponentController =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.205919, 1E-6);
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with LOCAL-MP2/TIGHT-PNO/Frozen core.
 */
TEST_F(MP2TaskTest, LocalMP2Restricted_TIGHT_FrozenCore) {
  SystemController__TEST_SUPPLY::cleanUp();
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.settings.nCoreOrbitals = 1;
  locTask.settings.useEnergyCutOff = false;
  locTask.run();
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(system, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  mp2Task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  mp2Task.settings.lcSettings.pnoThreshold = 1e-12;
  mp2Task.settings.lcSettings.orbitalToShellThreshold = 1e-5;
  mp2Task.settings.lcSettings.mullikenThreshold = 1e-5;
  mp2Task.settings.lcSettings.mullikenThreshold = 1e-5;
  mp2Task.settings.lcSettings.useFrozenCore = true;
  mp2Task.run();
  auto energyComponentController =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.20193476203093927, 1E-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests MP2Task.h/.cpp: Restricted energy test with LOCAL-MP2/Boughton--Pulay algorithm.
 */
TEST_F(MP2TaskTest, LocalMP2Restricted_TIGHT_BoughtonPulay) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  LocalizationTask locTask(system);
  locTask.run();
  MP2Task<Options::SCF_MODES::RESTRICTED> mp2Task(system, {});
  mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  mp2Task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  mp2Task.settings.lcSettings.pnoThreshold = 1e-12;
  mp2Task.settings.lcSettings.orbitalToShellThreshold = 1e-5;
  mp2Task.settings.lcSettings.mullikenThreshold = 1e-5;
  mp2Task.settings.lcSettings.useBPAlgorithm = true;
  mp2Task.settings.lcSettings.completenessThreshold = 0.0;
  mp2Task.run();
  auto energyComponentController =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  double MP2EnergyCorrection = energyComponentController->getEnergyComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
  EXPECT_NEAR(MP2EnergyCorrection, -0.205918, 1E-6);
}
} // namespace Serenity
