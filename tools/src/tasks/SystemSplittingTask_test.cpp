/**
 * @file SystemSplittingTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jan 9, 2020
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
#include "tasks/SystemSplittingTask.h"                //To be tested.
#include "data/ElectronicStructure.h"                 //DensityMatrix
#include "geometry/Geometry.h"                        //Default constructor.
#include "settings/Settings.h"                        //Settings.
#include "system/SystemController.h"                  //GetElectronicStructure.
#include "tasks/LocalizationTask.h"                   //Orbital localisation.
#include "tasks/ScfTask.h"                            //Clean electronic structures.
#include "tasks/SystemAdditionTask.h"                 //Supersystem construction.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class SystemSplittingTaskTest : public ::testing::Test {
 protected:
  SystemSplittingTaskTest() {
  }
  virtual ~SystemSplittingTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SystemSplittingTaskTest, restricted_mixedBasis) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  Settings settingsSys2 = sys2->getSettings();
  settingsSys2.basis.label = "DEF2-SVP";
  settingsSys2.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSys2);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.settings.addOccupiedOrbitals = false;
  add.run();

  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<RESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH;
  splitTask.run();

  EXPECT_EQ(0, sys1->getCharge());
  EXPECT_EQ(0, sys2->getCharge());

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(SystemSplittingTaskTest, restricted_populationThreshold) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto supersystem = *sys1 + *sys2;
  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<RESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD;
  // This is a extremely tiny value! This should select some orbitals of the second
  // water molecule for the first one too!
  splitTask.settings.orbitalThreshold = 0.001;
  splitTask.run();

  EXPECT_EQ(-2, sys1->getCharge());
  EXPECT_EQ(+2, sys2->getCharge());

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(SystemSplittingTaskTest, restricted_enforceCharges) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto supersystem = *sys1 + *sys2;
  sys1->setCharge(+2);
  sys2->setCharge(-2);
  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<RESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES;
  splitTask.run();

  EXPECT_EQ(+2, sys1->getCharge());
  EXPECT_EQ(4, sys1->getNOccupiedOrbitals<RESTRICTED>());
  EXPECT_EQ(-2, sys2->getCharge());
  EXPECT_EQ(6, sys2->getNOccupiedOrbitals<RESTRICTED>());

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(SystemSplittingTaskTest, restricted_enforceCharges2) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto supersystem = *sys1 + *sys2;
  sys1->setCharge(-2);
  sys2->setCharge(+2);
  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<RESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES;
  splitTask.run();

  EXPECT_EQ(-2, sys1->getCharge());
  EXPECT_EQ(6, sys1->getNOccupiedOrbitals<RESTRICTED>());
  EXPECT_EQ(+2, sys2->getCharge());
  EXPECT_EQ(4, sys2->getNOccupiedOrbitals<RESTRICTED>());

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(SystemSplittingTaskTest, restricted_SPADE) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto supersystem = *sys1 + *sys2;
  ScfTask<RESTRICTED> scf(supersystem);
  scf.run();

  SystemSplittingTask<RESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.settings.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::SPADE;
  splitTask.run();

  EXPECT_EQ(0, sys1->getCharge());
  EXPECT_EQ(5, sys1->getNOccupiedOrbitals<RESTRICTED>());
  EXPECT_EQ(0, sys2->getCharge());
  EXPECT_EQ(5, sys2->getNOccupiedOrbitals<RESTRICTED>());

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(SystemSplittingTaskTest, unrestricted_bestMatch) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::VERBOSE;
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto supersystem = *sys1 + *sys2;
  ScfTask<UNRESTRICTED> scf(supersystem);
  scf.run();

  LocalizationTask loc(supersystem);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.run();

  SystemSplittingTask<UNRESTRICTED> splitTask(supersystem, {sys1, sys2});
  splitTask.run();

  EXPECT_EQ(0, sys1->getCharge());
  EXPECT_EQ(0, sys2->getCharge());
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
