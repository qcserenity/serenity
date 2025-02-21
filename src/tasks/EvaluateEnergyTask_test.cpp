/**
 * @file EvaluateEnergyTask_test.cpp
 *
 * @date 6 Apr 2020
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
#include "tasks/EvaluateEnergyTask.h" //To be tested.
#include "data/ElectronicStructure.h" //GetEnergy.
#include "energies/EnergyContributions.h"
#include "settings/Settings.h"       //Get system with altered settings.
#include "system/SystemController.h" //GetElectronicStructure.
#include "tasks/FDETask.h"
#include "tasks/ScfTask.h"                            //Supersystem SCF.
#include "tasks/SystemAdditionTask.h"                 //Add fragments up.
#include "tasks/TDEmbeddingTask.h"                    //Run embedding calculation.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class EvaluateEnergyTaskTest : public ::testing::Test {
 protected:
  EvaluateEnergyTaskTest() {
  }
  virtual ~EvaluateEnergyTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(EvaluateEnergyTaskTest, restricted_compareSupersystemToExactEmbedding) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto superSettings = act->getSettings();
  auto super = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, superSettings);
  TDEmbeddingTask<RESTRICTED> tdTask(act, env);
  tdTask.run();
  double embeddingEnergy = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  SystemAdditionTask<RESTRICTED> addTask(super, {act, env});
  addTask.run();

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({super});
  evalEnergyTask.settings.XCfunctional = superSettings.dft.functional;
  evalEnergyTask.run();
  double superEnergy = super->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(embeddingEnergy, superEnergy, 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(super);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, restricted_sDFT) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);

  auto superSettings = act->getSettings();

  FDETask<RESTRICTED> fdeTask(act, {env});
  fdeTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  fdeTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  fdeTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  fdeTask.run();

  double energyBefore = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({act, env});
  evalEnergyTask.settings.XCfunctional = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  evalEnergyTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  evalEnergyTask.run();
  double energyAfter = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(energyBefore, energyAfter, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, unrestricted_sDFT) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE, true);

  auto superSettings = act->getSettings();

  FDETask<UNRESTRICTED> fdeTask(act, {env});
  fdeTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  fdeTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  fdeTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  fdeTask.run();

  double energyBefore = act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy();

  EvaluateEnergyTask<UNRESTRICTED> evalEnergyTask({act, env});
  evalEnergyTask.settings.XCfunctional = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  evalEnergyTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  evalEnergyTask.run();
  double energyAfter = act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy();
  EXPECT_NEAR(energyBefore, energyAfter, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, restricted_sDFTchangeXC) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);

  auto superSettings = act->getSettings();

  FDETask<RESTRICTED> fdeTask(act, {env});
  fdeTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  fdeTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  fdeTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  fdeTask.run();

  double energyBefore = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({act, env});
  evalEnergyTask.settings.useDifferentXCFunc = true;
  evalEnergyTask.settings.XCfunctional = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  evalEnergyTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  evalEnergyTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  evalEnergyTask.run();
  double energyAfter = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_TRUE((energyAfter < energyBefore - 1e-2) or (energyAfter > energyBefore + 1e-2));
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, restricted_sDFTchangeNaddXC) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);

  auto superSettings = act->getSettings();

  FDETask<RESTRICTED> fdeTask(act, {env});
  fdeTask.settings.embedding.naddXCFunc = superSettings.dft.functional;
  fdeTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  fdeTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  fdeTask.run();

  double energyBefore = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({act, env});
  evalEnergyTask.settings.XCfunctional = superSettings.dft.functional;
  evalEnergyTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  evalEnergyTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  evalEnergyTask.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  evalEnergyTask.run();
  double energyAfter = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_TRUE((energyAfter < energyBefore - 1e-4) or (energyAfter > energyBefore + 1e-4));
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, restricted_compareSystemToSCF_doubleHybrid_Local) {
  auto super = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  ScfTask<RESTRICTED> scfTask(super);
  scfTask.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  scfTask.run();
  double refEnergy =
      super->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({super});
  evalEnergyTask.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  evalEnergyTask.settings.lcSettings.diisStartResidual = 1e0;
  evalEnergyTask.settings.XCfunctional = (super->getSettings()).dft.functional;
  evalEnergyTask.run();
  double superEnergy =
      super->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
  EXPECT_NEAR(superEnergy, refEnergy, 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(EvaluateEnergyTaskTest, restricted_compareSystemToSCF_doubleHybrid_LocalVsRI) {
  auto super = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  ScfTask<RESTRICTED> scfTask(super);
  scfTask.settings.mp2Type = Options::MP2_TYPES::DF;
  scfTask.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  scfTask.run();
  double refEnergy = super->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();

  EvaluateEnergyTask<RESTRICTED> evalEnergyTask({super});
  evalEnergyTask.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  evalEnergyTask.settings.XCfunctional = (super->getSettings()).dft.functional;
  evalEnergyTask.run();
  double superEnergy = super->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(superEnergy, refEnergy, 5e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
