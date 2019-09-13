/**
 * @file CubeFileTask_test.cpp
 *
 * @date    Nov 24, 2015
 * @author: Jan Unsleber
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
#include "tasks/CubeFileTask.h"
#include "data/ElectronicStructure.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>


namespace Serenity {

/**
 * @class CubeFileTaskTest
 * @brief Sets everything up for the tests of CubeFileTask.h/.cpp .
 */
class CubeFileTaskTest : public ::testing::Test {
 protected:
  CubeFileTaskTest() {
  }

  virtual ~CubeFileTaskTest() = default;


  inline bool fileExists (const std::string& name) {
      std::ifstream f(name.c_str());
      if (f.good()) {
          f.close();
          return true;
      } else {
          f.close();
          return false;
      }
  }

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};


/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print density to file.
 */
TEST_F(CubeFileTaskTest, DensityRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.density =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_Density.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_Density.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print density to file with environment system.
 */
TEST_F(CubeFileTaskTest, DensityRestrictedEnvironmentSystem) {
  auto systemControllerA =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto systemControllerB =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemControllerA},{systemControllerB});
  task.settings.density =true;
  task.run();
  EXPECT_TRUE(fileExists((systemControllerA->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemControllerA->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print ELF to file.
 */
TEST_F(CubeFileTaskTest, ELFRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.elf =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_ELF.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_ELF.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print ELFTS to file.
 */
TEST_F(CubeFileTaskTest, ELFTSRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.elfts =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_ELF.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_ELF.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print SEDD to file.
 */
TEST_F(CubeFileTaskTest, SEDDRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.sedd =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_SEDD.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_SEDD.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print DORI to file.
 */
TEST_F(CubeFileTaskTest, DORIRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.dori =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_DORI.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_DORI.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print DORI to file.
 */
TEST_F(CubeFileTaskTest, SignedDensityRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.signedDensity =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_signedDensity.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_signedDensity.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print density to file.
 */
TEST_F(CubeFileTaskTest, DensityUnrestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto task =   CubeFileTask<UNRESTRICTED>({systemController},{});
  task.settings.density =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_alpha_Density.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_alpha_Density.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_beta_Density.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_beta_Density.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_TotalDensity.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_TotalDensity.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_SpinDensity.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_SpinDensity.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print MOs to file.
 */
TEST_F(CubeFileTaskTest, MOsRestricted) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.allOrbitals =true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO3.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print MOs to file.
 */
TEST_F(CubeFileTaskTest, MOsRestricted_parts) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.orbitals = {1};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO3.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print MOs to file.
 */
TEST_F(CubeFileTaskTest, MOsRestricted_nonExisiting) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.orbitals = {1,2,3};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_MO3.cube").c_str()));
}

/**
 * @test
 * @brief Tests CubeFileTask.h/.cpp: Print ESP to file.
 */
TEST_F(CubeFileTaskTest, ESP) {
  auto systemController =
    SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task =   CubeFileTask<RESTRICTED>({systemController},{});
  task.settings.electrostaticPot = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSettings().path+"TestSystem_H2_MINBAS_ESP.cube").c_str()));
  EXPECT_EQ(0,std::remove((systemController->getSettings().path+"TestSystem_H2_MINBAS_ESP.cube").c_str()));
}



} /* namespace Serenity */
