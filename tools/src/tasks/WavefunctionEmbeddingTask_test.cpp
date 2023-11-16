/**
 * @file WavefunctionEmbeddingTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 29, 2020
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "tasks/WavefunctionEmbeddingTask.h" //To be tested.
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h" //Default geometry constructor.
#include "settings/Settings.h" //Supersystem construction.
#include "system/SystemController.h"
#include "tasks/CoupledClusterTask.h"                 //Supersystem CCSD.
#include "tasks/LocalizationTask.h"                   //Orbital localization.
#include "tasks/ScfTask.h"                            //Supersystem scf.
#include "tasks/SystemAdditionTask.h"                 //Supersystem construction.
#include "tasks/SystemSplittingTask.h"                //System partitioning.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class WavefunctionEmbeddingTaskTest : public ::testing::Test {
 protected:
  WavefunctionEmbeddingTaskTest() = default;

  virtual ~WavefunctionEmbeddingTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
  void cleanUp(std::shared_ptr<SystemController> system) {
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  }
};

TEST_F(WavefunctionEmbeddingTaskTest, waterDimer_supersystemSplitToFragments) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  Settings settings = A->getSettings();
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settings);
  SystemAdditionTask<RESTRICTED> add(supersystem, {A, B});
  add.settings.addOccupiedOrbitals = false;
  add.run();

  WavefunctionEmbeddingTask wfEmb(supersystem, {A, B});
  wfEmb.settings.normThreshold = 1e-9;
  wfEmb.settings.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::NORMAL;
  wfEmb.settings.lcSettings[0].method = Options::PNO_METHOD::DLPNO_CCSD;
  wfEmb.settings.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::NORMAL;
  wfEmb.settings.lcSettings[1].method = Options::PNO_METHOD::DLPNO_CCSD;
  wfEmb.settings.fromFragments = false;
  wfEmb.settings.loc.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  wfEmb.settings.loc.splitValenceAndCore = true;
  wfEmb.settings.split.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH;
  wfEmb.run();

  CoupledClusterTask ccTask(supersystem, {});
  ccTask.settings.level = Options::CC_LEVEL::DLPNO_CCSD;
  ccTask.settings.normThreshold = 1e-9;
  ccTask.run();

  auto fragmentEnergies = wfEmb.getFragmentEnergies();
  auto interactionEnergies = wfEmb.getFragmentWiseInteractionEnergy();

  EXPECT_NEAR(fragmentEnergies.sum() + interactionEnergies(0),
              supersystem->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION), 1e-9);

  cleanUp(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(WavefunctionEmbeddingTaskTest, waterDimer_supersystemFromFragments) {
  auto A = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  auto B = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  Settings settings = A->getSettings();
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settings);
  SystemAdditionTask<RESTRICTED> add(supersystem, {A, B});
  add.settings.addOccupiedOrbitals = false;
  add.run();

  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.splitValenceAndCore = true;
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<RESTRICTED> split(supersystem, {A, B});
  split.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH;
  split.run();

  WavefunctionEmbeddingTask wfEmb(supersystem, {A, B});
  wfEmb.settings.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::NORMAL;
  wfEmb.settings.lcSettings[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  wfEmb.settings.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::NORMAL;
  wfEmb.settings.lcSettings[1].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  wfEmb.settings.fromFragments = true;
  wfEmb.settings.normThreshold = 1e-9;
  wfEmb.run();

  CoupledClusterTask ccTask(supersystem, {});
  ccTask.settings.level = Options::CC_LEVEL::DLPNO_CCSD_T0;
  ccTask.settings.normThreshold = 1e-9;
  ccTask.run();

  auto fragmentEnergies = wfEmb.getFragmentEnergies();
  auto interactionEnergies = wfEmb.getFragmentWiseInteractionEnergy();

  EXPECT_NEAR(fragmentEnergies.sum() + interactionEnergies(0),
              supersystem->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION) +
                  supersystem->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION),
              1e-9);

  cleanUp(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
