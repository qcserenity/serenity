/**
 * @file PCMInteractionEnergyTask_test.cpp
 *
 * @date   Nov 12, 2020
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
#include "tasks/PCMInteractionEnergyTask.h" //To be tested.
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class PCMInteractionEnergyTaskTest : public ::testing::Test {
 protected:
  PCMInteractionEnergyTaskTest() {
  }

  virtual ~PCMInteractionEnergyTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(PCMInteractionEnergyTaskTest, cpcm_water) {
  auto water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);

  PCMInteractionEnergyTask task(water);
  task.settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  task.settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  task.run();

  auto energy = water->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY);
  EXPECT_NEAR(-0.00946404, energy, 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
