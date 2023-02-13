/**
 * @file FDETask_test.cpp
 *
 * @date Mar 14, 2017
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
#include "tasks/FDETask.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class FDETask
 * @brief Sets everything up for the tests of FDETask.h/.cpp .
 */
class FDETaskTest : public ::testing::Test {
 protected:
  FDETaskTest() {
  }

  virtual ~FDETaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FDETaskTest, restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FDETask<Options::SCF_MODES::RESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-1.8129554924421183, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FDETaskTest, restricted_cut_grid) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FDETask<Options::SCF_MODES::RESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.gridCutOff = 5.0;
  task.run();
  EXPECT_NEAR(-1.8129555447384618, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy.
 */
TEST_F(FDETaskTest, unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  auto task = FDETask<Options::SCF_MODES::UNRESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-1.8081854183941235, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FDETaskTest, restricted_unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  auto task = FDETask<Options::SCF_MODES::RESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-1.8081854183941235, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy.
 */
TEST_F(FDETaskTest, unrestricted_restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(RESTRICTED);
  auto task = FDETask<Options::SCF_MODES::UNRESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-1.8129554921695687, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test PCM embedding.
 */
TEST_F(FDETaskTest, restrictedPCM) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FDETask<Options::SCF_MODES::RESTRICTED>(act, {env});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.pcm.use = true;
  task.settings.embedding.pcm.solvent = Options::PCM_SOLVENTS::N_HEPTANE;
  task.settings.embedding.pcm.solverType = Options::PCM_SOLVER_TYPES::IEFPCM;
  task.run();
  EXPECT_NEAR(-1.8131957184036818, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test the unrelaxed MP2 density correction.
 */
TEST_F(FDETaskTest, densityCorrectedInteraction) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  auto task = FDETask<Options::SCF_MODES::RESTRICTED>(act, {env});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.calculateUnrelaxedMP2Density = true;
  task.run();
  double interaction_energy =
      act->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY);
  EXPECT_NEAR(-0.50363547208851855, interaction_energy, 1e-7);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
