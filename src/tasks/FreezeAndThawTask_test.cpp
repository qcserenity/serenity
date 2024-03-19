/**
 * @file FreezeAndThawTask_test.cpp
 *
 * @date Mar 16, 2017
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
#include "tasks/FreezeAndThawTask.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "io/IOOptions.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class FDETask
 * @brief Sets everything up for the tests of ExportGridTask.h/.cpp .
 */
class FreezeAndThawTaskTest : public ::testing::Test {
 protected:
  FreezeAndThawTaskTest() {
  }

  virtual ~FreezeAndThawTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  double getSupersystemEnergy() {
    if (_supersystemEnergy < 1e-10) {
      auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
      auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
      auto supersystem = *act + *env;
      ScfTask<Options::SCF_MODES::RESTRICTED> task(supersystem);
      task.run();
      _supersystemEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
      // clean up
      SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
    } // if
    return _supersystemEnergy;
  }

 private:
  double _supersystemEnergy = 0.0;
};

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FreezeAndThawTaskTest, restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});

  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.run();
  EXPECT_NEAR(-1.7189271379482964, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FreezeAndThawTaskTest, restrictedPassive) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act}, {env});

  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.run();
  EXPECT_NEAR(-1.7139541221231887, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FreezeAndThawTaskTest, restricted_cut_grid) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});

  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.gridCutOff = 5.0;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.run();
  EXPECT_NEAR(-1.7189270899267217, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy.
 */
TEST_F(FreezeAndThawTaskTest, unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({act, env}, {});
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.run();
  EXPECT_NEAR(-1.812956512024005, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with Hoffmann's projection technique.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the Huzinaga equation and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Huzinaga_DIIS) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  iOOptions.printDebugInfos = true;
  task.run();
  iOOptions.printDebugInfos = false;

  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the Huzinaga equation and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_HYBRID) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  task.run();
  auto supersystem = *act + *env;
  ScfTask<Options::SCF_MODES::RESTRICTED> SuperScf(supersystem);
  SuperScf.run();
  double superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  // clean up
  auto supersystemName = supersystem->getSystemName();
  EXPECT_NEAR(superEnergy, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem->getSystemPath() + "/", supersystemName);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with Hoffmann's projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_DIIS) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy with Hoffmann's projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_DIIS_UNRESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({act, env}, {});
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the level-shift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Levelshift_DIIS) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with HF-in-DFT, Levelshift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_WaterDimer_HFinDFT_Levelshift_DIIS) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(-152.4201802840, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy after TD calculation.
 */
TEST_F(FreezeAndThawTaskTest, TDPlusFaT) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  // Perform TD calculation
  TDEmbeddingTask<SPIN> tdTask(act, env);
  tdTask.run();
  // Perform FaT calculation on top
  auto task = FreezeAndThawTask<SPIN>({act, env}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.run();
  // Different results to the scratch calculation (restricted) due to different settings of the systems on disk.
  EXPECT_NEAR(getSupersystemEnergy(), act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  // Clean up
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with HF-in-DFT, Levelshift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, WaterDimer_HFinDFT_PassiveSolvated) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
  task.settings.calculateSolvationEnergy = true;
  task.settings.finalEnergyEvaluation = false;
  task.run();

  auto eCont = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  EXPECT_NEAR(-76.012318767499409, eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with HF-in-DFT, Levelshift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, WaterDimer_DFTinDFT_PassiveSolvated) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
  task.settings.calculateSolvationEnergy = true;
  task.settings.finalEnergyEvaluation = false;
  task.run();

  auto eCont = act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  EXPECT_NEAR(-76.411149223007001,
              eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Tests energy for multiscale approx in exact embedding (Levelshift/ nadd-kin).
 */
TEST_F(FreezeAndThawTaskTest, rMultiscaleLevel) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act1_act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86, true);
  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86, true);
  auto env3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86, true);
  act1_act2->setSCFMode(SPIN);
  act1->setSCFMode(SPIN);
  act2->setSCFMode(SPIN);
  env3->setSCFMode(SPIN);
  auto task = FreezeAndThawTask<SPIN>({act1_act2, env3}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  // Perform FaT with comibed 1+2 and 3 calculation
  task.run();
  // Perform Localization
  auto localization = LocalizationTask(act1_act2);
  localization.run();
  // Splitting
  auto splitting = SystemSplittingTask<SPIN>(act1_act2, {act1, act2});
  splitting.run();
  // Run additional FaT with 1 cycle to relax exactly treated systems
  // auto task2 = FDETask<SPIN>(act1,{act2,env3});
  // task2.settings.calculateEnvironmentEnergy = true;
  auto task2 = FreezeAndThawTask<SPIN>({act1, act2}, {env3});
  task2.settings.maxCycles = 2;
  task2.settings.convThresh = 0.0;
  task2.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task2.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                             CompositeFunctionals::XCFUNCTIONALS::BP86};
  task2.run();
  EXPECT_NEAR(act1_act2->getElectronicStructure<SPIN>()->getEnergy(), act2->getElectronicStructure<SPIN>()->getEnergy(), 2e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Tests energy for multiscale approx in exact embedding (Huzinaga/ nadd-kin; unrestricted).
 */
TEST_F(FreezeAndThawTaskTest, uMultiscaleHuzinaga) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;
  auto act1_act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86, true);
  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86, true);
  auto env3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86, true);
  act1_act2->setSCFMode(SPIN);
  act1->setSCFMode(SPIN);
  act2->setSCFMode(SPIN);
  env3->setSCFMode(SPIN);
  auto task = FreezeAndThawTask<SPIN>({act1_act2, env3}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  // Perform FaT with comibed 1+2 and 3 calculation
  task.run();
  // Perform Localization
  auto localization = LocalizationTask(act1_act2);
  localization.run();
  // Splitting
  auto splitting = SystemSplittingTask<SPIN>(act1_act2, {act1, act2});
  splitting.run();
  // Run additional FaT with 1 cycle to relax exactly treated systems
  auto task2 = FreezeAndThawTask<SPIN>({act1, act2}, {env3});
  task2.settings.maxCycles = 2;
  task2.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::HUZINAGA, Options::KIN_EMBEDDING_MODES::HUZINAGA,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task2.settings.convThresh = 0.0;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task2.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                             CompositeFunctionals::XCFUNCTIONALS::BP86};
  task2.run();

  EXPECT_NEAR(act1_act2->getElectronicStructure<SPIN>()->getEnergy(), act2->getElectronicStructure<SPIN>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}
} /*namespace Serenity*/
