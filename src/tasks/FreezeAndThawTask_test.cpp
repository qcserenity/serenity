/**
 * @file FreezeAndThawTask_test.cpp
 *
 * @date Mar 16, 2017
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

/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "tasks/FreezeAndThawTask.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "tasks/ScfTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "io/IOOptions.h"
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
      auto act =
        SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
      auto env =
          SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
      auto supersystem = *act + *env;
      ScfTask<Options::SCF_MODES::RESTRICTED>task (supersystem);
      task.run();
      _supersystemEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
    // clean up
    std::remove((supersystem->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
    std::remove((supersystem->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
    std::remove((supersystem->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
    std::remove((supersystem->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
    std::remove((supersystem->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
    std::remove((supersystem->getSettings().path).c_str());
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
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE,true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});

  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 0;
  task.run();
  EXPECT_NEAR(-1.7189271379482964,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();



}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy.
 */
TEST_F(FreezeAndThawTaskTest, restricted_cut_grid) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE,true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});

  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.gridCutOff = 5.0;
  task.settings.printLevel = 0;
  task.run();
  EXPECT_NEAR(-1.7189270899267217,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();



}


/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy.
 */
TEST_F(FreezeAndThawTaskTest, unrestricted) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({act,env},{});
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 1;
  task.run();
  EXPECT_NEAR(-1.812956512024005,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({act,env},{});
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.settings.makeSuperSystemBasis = true;
  task.run();
  EXPECT_NEAR(-1.8190429139113393,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with Hoffmann's projection technique.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the Huzinaga equation and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Huzinaga_DIIS) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  iOOptions.printDebugInfos = true;
  task.run();
  iOOptions.printDebugInfos = false;

  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the Huzinaga equation and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_HYBRID) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID,true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID,true);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.basisExtThresh = 0.0;
  task.settings.extendBasis = true;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
  task.settings.printLevel = 2;
  task.run();
  auto supersystem = *act + *env;
  ScfTask<Options::SCF_MODES::RESTRICTED>SuperScf (supersystem);
  SuperScf.run();
  double superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  // clean up
  auto supersystemName=supersystem->getSystemName();
  std::remove((supersystem->getSettings().path+supersystemName+".settings").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".xyz").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".energies.res.h5").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".orbs.res.h5").c_str());
  std::remove((supersystem->getSettings().path+supersystemName+".dmat.res.h5").c_str());
  std::remove((supersystem->getSettings().path).c_str());

  EXPECT_NEAR(superEnergy,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}


/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with Hoffmann's projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_DIIS) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test unrestricted energy with Hoffmann's projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Hoffmann_DIIS_UNRESTRICTED) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({act,env},{});
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with the level-shift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_Levelshift_DIIS) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.run();

  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy with HF-in-DFT, Levelshift projection technique and DIIS.
 */
TEST_F(FreezeAndThawTaskTest, supersystembasis_WaterDimer_HFinDFT_Levelshift_DIIS) {
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);
  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act,env},{});
  act->setSCFMode(RESTRICTED);
  env->setSCFMode(RESTRICTED);
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.useConvAcceleration = true;
  task.settings.diisStart = 0.1;
  task.settings.extendBasis = true;
  task.settings.basisExtThresh = 0.0;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 2;
  task.run();

  EXPECT_NEAR(-152.4201802840,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests FDETask.h/.cpp: Test restricted energy after TD calculation.
 */
TEST_F(FreezeAndThawTaskTest, TDPlusFaT) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  // Perform TD calculation
  TDEmbeddingTask<SPIN> tdTask(act,env);
  tdTask.run();
  // Perform FaT calculation on top
  auto task = FreezeAndThawTask<SPIN>({act,env},{});
  task.settings.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.printLevel = 0;
  task.run();
  // Different results to the scratch calculation (restricted) due to different settings of the systems on disk.
  EXPECT_NEAR(getSupersystemEnergy(),act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);
  // Clean up
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
