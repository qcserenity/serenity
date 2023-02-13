/**
 * @file GeneralizedDOSTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 24, 2020
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
#include "tasks/GeneralizedDOSTask.h"                               //To be tested.
#include "analysis/directOrbitalSelection/DirectOrbitalSelection.h" //Orbital groups definition
#include "geometry/Geometry.h"                                      //Check geometries after selection.
#include "settings/Settings.h"                                      //Settings definition.
#include "system/SystemController.h"                                //Access to test systems.
#include "tasks/LocalizationTask.h"                                 //Orbital localization.
#include "tasks/PlotTask.h"
#include "tasks/WavefunctionEmbeddingTask.h"          //Coupled cluster calculation.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class GeneralizedDOSTaskTest : public ::testing::Test {
 protected:
  GeneralizedDOSTaskTest() = default;

  virtual ~GeneralizedDOSTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(GeneralizedDOSTaskTest, restricted_EthaneIAOShell) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);

  Settings fragSettings = a->getSettings();
  fragSettings.name = "a1";
  auto a1 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "a2";
  auto a2 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "a3";
  auto a3 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b1";
  auto b1 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b2";
  auto b2 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b3";
  auto b3 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;

  LocalizationTask loc(a);
  loc.settings.splitValenceAndCore = true;
  loc.settings.localizeVirtuals = true;
  loc.run();

  LocalizationTask loc2(b);
  loc2.settings = loc.settings;
  loc2.run();

  GeneralizedDOSTask<SPIN> gdos({a, b}, {a1, a2, a3, b1, b2, b3});
  gdos.settings.similarityLocThreshold = {1e-1, 5e-2};
  gdos.settings.similarityKinEnergyThreshold = {1e-1, 5e-2};
  gdos.settings.usePiBias = true;
  gdos.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell;
  gdos.settings.mapVirtuals = true;
  gdos.settings.writeGroupsToFile = true;
  gdos.run();

  EXPECT_EQ(1, a1->getNOccupiedOrbitals<SPIN>());
  EXPECT_EQ(6, a2->getNOccupiedOrbitals<SPIN>());
  EXPECT_EQ(2, a3->getNOccupiedOrbitals<SPIN>());
  EXPECT_EQ(1, b1->getNOccupiedOrbitals<SPIN>());
  EXPECT_EQ(6, b2->getNOccupiedOrbitals<SPIN>());
  EXPECT_EQ(2, b3->getNOccupiedOrbitals<SPIN>());

  const auto& groups = gdos.getOrbitalGroups();
  for (const auto& group : groups) {
    EXPECT_EQ(1, group->getNOrbitals());
  }

  GeneralizedDOSTask<SPIN> gdos2({a, b}, {a1, a2, a3, b1, b2, b3});
  gdos2.settings = gdos.settings;
  gdos2.settings.bestMatchMapping = true;
  gdos2.settings.scoreStart = 5e-1;
  gdos2.settings.scoreEnd = 5e-2;
  gdos2.run();
  const auto& groups2 = gdos2.getOrbitalGroups();
  for (const auto& group : groups2) {
    if (group->getNOrbitals()) {
      EXPECT_EQ(1, group->getNOrbitals());
    }
  }

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;

  std::remove("orbital_groups.dat");

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a3);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b3);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(GeneralizedDOSTaskTest, restricted_gdos_wfemb_run) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86);

  // Switch to HF
  Settings fragSettings = a->getSettings();
  fragSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  fragSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, fragSettings);
  fragSettings.name = "b";
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, fragSettings);
  // Build fragments.
  fragSettings.name = "a1";
  auto a1 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "a2";
  auto a2 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "a3";
  auto a3 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b1";
  auto b1 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b2";
  auto b2 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  fragSettings.name = "b3";
  auto b3 = std::make_shared<SystemController>(std::make_shared<Geometry>(), fragSettings);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;

  LocalizationTask aLoc(a, {});
  aLoc.settings.splitValenceAndCore = true;
  aLoc.run();
  LocalizationTask bLoc(b, {});
  bLoc.settings.splitValenceAndCore = true;
  bLoc.run();

  GeneralizedDOSTask<SPIN> gdos({a, b}, {a1, a2, a3, b1, b2, b3});
  gdos.settings.similarityLocThreshold = {1e-1, 5e-2};
  gdos.settings.similarityKinEnergyThreshold = {1e-1, 5e-2};
  gdos.settings.usePiBias = true;
  gdos.settings.populationAlgorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell;
  gdos.run();

  WavefunctionEmbeddingTask wfemb(a, {a1, a2, a3});
  wfemb.settings.lcSettings[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  wfemb.settings.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  wfemb.settings.lcSettings[1].method = Options::PNO_METHOD::DLPNO_CCSD;
  wfemb.settings.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::NORMAL;
  wfemb.settings.lcSettings[2].method = Options::PNO_METHOD::DLPNO_CCSD;
  wfemb.settings.lcSettings[2].pnoSettings = Options::PNO_SETTINGS::LOOSE;
  wfemb.settings.fromFragments = true;
  wfemb.settings.fullDecomposition = false;
  wfemb.run();

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a3);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b3);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
