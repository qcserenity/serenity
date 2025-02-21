/**
 * @file TopDownStaticEmbeddingTask_test.cpp
 *
 * @date Mar 21, 2024
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
#include "tasks/TopDownStaticEmbeddingTask.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class TopDownStaticEmbeddingTaskTest : public ::testing::Test {
 protected:
  TopDownStaticEmbeddingTaskTest() {
  }

  virtual ~TopDownStaticEmbeddingTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests TopDownStaticEmbeddingTask.h/.cpp: Test restricted energy for multiple subsystems with exact DFT-in-DFT
 * embedding.
 */
TEST_F(TopDownStaticEmbeddingTaskTest, multipleSubsystemsDFT) {
  auto actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, true);
  auto actB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B, true);
  auto actC = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C, true);
  auto envD = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D, true);
  auto envE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_E, true);
  auto envF = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_F, true);

  TopDownStaticEmbeddingTask<RESTRICTED> task({actA, actB, actC}, {envD, envE, envF});
  task.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.lcSettings.embeddingSettings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  task.run();

  auto totalEnergy = task.getFinalEnergy();
  auto supersystemReferenceEnergy = -457.7165269246;
  EXPECT_NEAR(totalEnergy, supersystemReferenceEnergy, 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envD->getSystemPath() + "TMP_Supersystem/", "TMP_Supersystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TopDownStaticEmbeddingTask.h/.cpp: Test restricted energy for multiple subsystems with MP2-in-DFT
 * embedding. The reference is only taken from a "working" implementation.
 */
TEST_F(TopDownStaticEmbeddingTaskTest, multipleSubsystemsMP2inDFT) {
  auto actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, true);
  auto actB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B, true);
  auto actC = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C, true);
  actA->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::HF);
  actB->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::HF);
  actC->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::HF);
  auto envD = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D, true);
  auto envE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_E, true);
  auto envF = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_F, true);

  TopDownStaticEmbeddingTask<RESTRICTED> task({actA, actB, actC}, {envD, envE, envF});
  task.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.lcSettings.embeddingSettings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  task.settings.lcSettings.maximumMemoryRatio = 0.0;
  task.settings.lcSettings.ignoreMemoryConstraints = true;
  task.run();

  auto totalEnergy = task.getFinalEnergy();
  auto supersystemReferenceEnergy = -457.37706342397757;
  EXPECT_NEAR(totalEnergy, supersystemReferenceEnergy, 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envD->getSystemPath() + "TMP_Supersystem/", "TMP_Supersystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(TopDownStaticEmbeddingTaskTest, multipleSubsystemsDFTinDFTExternalCharges) {
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/TestSystem_Water_Hexamer_A/external_charges.pc";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }

  auto actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, true);
  auto actB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B, true);
  auto actC = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C, true);
  auto envD = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D, true);
  std::vector<std::shared_ptr<SystemController>> allSubsystems = {actA, actB, actC, envD};
  std::vector<Settings> allSettings;
  for (auto& subsystem : allSubsystems) {
    Settings settings = subsystem->getSettings();
    settings.extCharges.externalChargesFile = pathToTestsResources;
    allSettings.push_back(settings);
  }
  actA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_A, allSettings[0]);
  actB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_B, allSettings[1]);
  actC = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_C, allSettings[2]);
  envD = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Hexamer_Monomer_D, allSettings[3]);

  TopDownStaticEmbeddingTask<RESTRICTED> task({actA, actB, actC}, {envD});
  task.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.lcSettings.embeddingSettings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  task.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  task.run();

  auto totalEnergy = task.getFinalEnergy();
  auto supersystemReferenceEnergy = -305.1617561304;
  EXPECT_NEAR(totalEnergy, supersystemReferenceEnergy, 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(envD->getSystemPath() + "TMP_Supersystem/", "TMP_Supersystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */