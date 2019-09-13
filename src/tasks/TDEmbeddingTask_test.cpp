/**
 * @file TDEmbeddingTask_test.cpp
 *
 * @date Mar 14, 2017
 * @author Jan Unsleber
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

#include "data/ElectronicStructure.h"
#include "tasks/TDEmbeddingTask.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
  /* Include Serenity Internal Headers */
/**
 * @class TDEmbeddingTaskTest
 * @brief Sets everything up for the tests of TDEmbeddingTask.h/.cpp .
 */
class TDEmbeddingTaskTest : public ::testing::Test {
 protected:
  TDEmbeddingTaskTest() {
  }
  virtual ~TDEmbeddingTaskTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy.
 */
TEST_F(TDEmbeddingTaskTest, restricted) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.run();
  EXPECT_NEAR(-1.8155882075161134,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy.
 */
TEST_F(TDEmbeddingTaskTest, restricted_read_supersystem) {
    std::string pathToTestsResources;
  if(const char* env_p = std::getenv("SERENITY_RESOURCES")){
    pathToTestsResources = (std::string)env_p + "testresources/";
  }else{
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set."<< std::endl;
    assert(false);

  }
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.load = pathToTestsResources;
  task.settings.name = "TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE";
  task.run();
  EXPECT_NEAR(-1.8155882075161134,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE/TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE.settings");
  std::remove("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE/TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE.xyz");
  std::remove("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE/TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE.dmat.res.h5");
  std::remove("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act,env);
  task.run();
  EXPECT_NEAR(-1.8155882296534891,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy when relaxing to precalculated environment.
 */
TEST_F(TDEmbeddingTaskTest, useEnvSys) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.enforceCharges = true;
  task.run();



  auto task2 = TDEmbeddingTask<RESTRICTED>(env,act);
  task2.settings.useEnvSys = true;
  task2.settings.enforceCharges = true;
  task2.run();

  std::remove("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE/TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE.orbs.res.h5");
  EXPECT_NEAR(-1.815588208332326,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE.settings").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE.xyz").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE.energies.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE.orbs.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE.dmat.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE/out").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE+TestSystem_H2_6_31Gs_ACTIVE_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy with Hoffmann's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_hoffmann) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.run();
  EXPECT_NEAR(-1.8155881852554803,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy with Huzinaga's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_huzinaga) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-1.8155881852542037,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_huzinaga) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-1.8155881853192255,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with a fermi-shifted Huzinaga projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_fermi_huzinaga) {
  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-306.6721662019,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);
  EXPECT_NEAR(task.settings.embedding.fermiShift,3.430175862028e-02,1e-7);

  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with a shifted Huzinaga projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_shifted_huzinaga) {
  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.useFermiLevel = false;
  task.run();
  EXPECT_NEAR(-306.6721662019,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);
  EXPECT_NEAR(task.settings.embedding.fermiShift,1.0,1e-7);

  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga's projection operator for hybrid functionals.
 */
TEST_F(TDEmbeddingTaskTest, restricted_huzinaga_hybrid) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID,true);

  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  auto supersystem = *act+*env;
  ScfTask<Options::SCF_MODES::RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.run();
  auto superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(superEnergy,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_HYBRID+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_HYBRID").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT,true);
  // HF environment
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.netThreshold = 1e-4;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();
  EXPECT_NEAR(-152.42058006269133,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs").c_str());
	SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation and subsequent restart.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_restart) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT,true);
  // HF environment
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.netThreshold = 1e-4;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load=".";
  loadSettings.path="new";

  auto newSys=std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(),newSys->getBasisController()->getNBasisFunctions());

  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs").c_str());

  std::remove((newSys->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT.settings").c_str());
  std::remove((newSys->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT.xyz").c_str());
  std::remove((newSys->getSettings().path).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_restart_basisChange) {
  auto env =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS,true);
  // HF environment
  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.netThreshold = 1e-3;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load=".";
  loadSettings.path="new";

  auto newSys=std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(),newSys->getBasisController()->getNBasisFunctions());
  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(),(unsigned int)18);

  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS").c_str());

  std::remove((newSys->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT.settings").c_str());
  std::remove((newSys->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT.xyz").c_str());
  std::remove((newSys->getSettings().path).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_levelshift_restart_basisChange) {
  auto env =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS,true);
  // HF environment
  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load=".";
  loadSettings.path="new";

  auto newSys=std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(),newSys->getBasisController()->getNBasisFunctions());
  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(),(unsigned int)36);

  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS").c_str());

  std::remove((newSys->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT.settings").c_str());
  std::remove((newSys->getSettings().path+"TestSystem_WaterMonTwo_6_31Gs_DFT.xyz").c_str());
  std::remove((newSys->getSettings().path).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + prim. net population.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_primNetPop_levelshift) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.truncationFactor = 0.5;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();
  EXPECT_NEAR(-1.8133391420948282,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_huzinaga) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.netThreshold = 8e-2;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-1.8189844336164807,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population + non ortho. distant + kin. energy functional.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_nonortho_distant_kinetic) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT,true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.netThreshold = 1e-3;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.longRangeNaddKinFunc = Options::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.borderAtomThreshold = 0.8;
  task.settings.embedding.basisFunctionRatio = 0.4;
  task.run();
  EXPECT_NEAR(-152.8156966216045,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs_DFT").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy for the Carter potential reconstruction.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_Carter) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::RECONSTRUCTION;
  task.settings.embedding.carterCycles = 5000;
  task.run();
  EXPECT_NEAR(-1.8155881842923347,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population + non ortho. distant + kin. energy functional.
 */
TEST_F(TDEmbeddingTaskTest, restricted_WYLB) {

  auto act =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::RECONSTRUCTION;
  task.settings.embedding.lbCycles = 5000;
  task.settings.noSupRec = false;
  task.run();
  EXPECT_NEAR(-1.8062078106336934,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-5);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(TDEmbeddingTaskTest, embeddingMode_NONE) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  task.run();
  EXPECT_NEAR(-2.4352547558648183,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(TDEmbeddingTaskTest, embeddingMode_NADDFUNC) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.run();
  EXPECT_NEAR(-1.7947783381675915,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-7);

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests the TDEmbedding with ECPs.
 */
TEST_F(TDEmbeddingTaskTest, embeddingMode_Huzinaga_ECPs) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  // There is (currently) no MINAO basis for I.
  task.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE;
  task.run();
  // Perform supersystem calculations
  auto supersystem = *act+*env;
  ScfTask<Options::SCF_MODES::RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.run();
  auto superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(superEnergy,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-4);
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE").c_str());
  // Clean the remains of the supersystem calculation up.
  auto supersystemName = supersystem->getSystemName();
  std::remove((supersystem->getSettings().path+supersystemName+".settings").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".xyz").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".orbs.res.h5").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".dmat.res.h5").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".energies.res.h5").c_str());
  std::remove((supersystem->getSettings().path).c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the TDEmbedding with ECPs.
 */
TEST_F(TDEmbeddingTaskTest, embeddingMode_HuzinagaTruncated_ECPs) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  // There is (currently) no MINAO basis for I.
  task.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.netThreshold = 0.0;
  task.run();
  EXPECT_NEAR(-893.6684430205,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-4);
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
