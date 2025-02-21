/**
 * @file ImportCavityTask_test.cpp
 *
 * @date   Nov 25, 2024
 * @author Lukas Paetow
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
#include "tasks/ImportCavityTask.h" //To be tested.
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/ExportCavityTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ImportCavityTaskTest : public ::testing::Test {
 protected:
  ImportCavityTaskTest() {
  }

  virtual ~ImportCavityTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ImportCavityTaskTest, h2) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  Settings settings = systemController->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  settings.pcm.cavityFormation = true;
  settings.pcm.saveCharges = true;
  systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto es = systemController->getElectronicStructure<RESTRICTED>();

  ExportCavityTask task(systemController);
  task.run();

  // importing
  systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);

  ImportCavityTask task2(systemController);
  task2.settings.cavityPath = (systemController->getSystemPath());
  task2.settings.vdwcavityPath = (systemController->getSystemPath());

  task2.run();

  auto energy = systemController->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY);
  EXPECT_NEAR(0.0053495458, energy, 1e-7);

  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/CavityData.h5").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/VDWCavityData.h5").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/PCMChargesData.h5").c_str()));
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
