/**
 * @file PlotTask_test.cpp
 *
 * @date    Nov 24, 2015
 * @author: Jan Unsleber
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
#include "tasks/PlotTask.h"
#include "data/ElectronicStructure.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/Analysis/NROCalculator.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <string>

namespace Serenity {

/**
 * @class PlotTaskTest
 * @brief Sets everything up for the tests of PlotTask.h/.cpp .
 */
class PlotTaskTest : public ::testing::Test {
 protected:
  PlotTaskTest() {
  }

  virtual ~PlotTaskTest() = default;

  inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    if (f.good()) {
      f.close();
      return true;
    }
    else {
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
 * @brief Tests PlotTask.h/.cpp: Print density to file.
 */
TEST_F(PlotTaskTest, DensityRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.density = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density.cube").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print density to file.
 */
TEST_F(PlotTaskTest, DensityRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.density = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print density to heatmap.
 */
TEST_F(PlotTaskTest, DensityRestrictedHeatmap) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.density = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.settings.xyHeatmap = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density_XYPLANE.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density_XYPLANE.dat").c_str()));
  EXPECT_TRUE(fileExists(
      (systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density_MOLECULE_ROTATED_TO_XYPLANE.xyz").c_str()));
  EXPECT_EQ(0, std::remove(
                   (systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density_MOLECULE_ROTATED_TO_XYPLANE.xyz").c_str()));
  std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_Density.dat").c_str());
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print density to file with environment system.
 */
TEST_F(PlotTaskTest, DensityRestrictedEnvironmentSystemCube) {
  auto systemControllerA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto systemControllerB =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemControllerA}, {systemControllerB});
  task.settings.density = true;
  task.run();
  EXPECT_TRUE(fileExists((systemControllerA->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemControllerA->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.cube").c_str()));
}

TEST_F(PlotTaskTest, DensityRestrictedEnvironmentSystemPlane) {
  auto systemControllerA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto systemControllerB =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemControllerA}, {systemControllerB});
  task.settings.density = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemControllerA->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemControllerA->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE_Density.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print ELF to file.
 */
TEST_F(PlotTaskTest, ELFRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.elf = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.cube").c_str()));
}
TEST_F(PlotTaskTest, ELFRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.elf = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print ELFTS to file.
 */
TEST_F(PlotTaskTest, ELFTSRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.elfts = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.cube").c_str()));
}
TEST_F(PlotTaskTest, ELFTSRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.elfts = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ELF.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print SEDD to file.
 */
TEST_F(PlotTaskTest, SEDDRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.sedd = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SEDD.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SEDD.cube").c_str()));
}
TEST_F(PlotTaskTest, SEDDRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.sedd = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SEDD.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SEDD.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print DORI to file.
 */
TEST_F(PlotTaskTest, DORIRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.dori = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_DORI.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_DORI.cube").c_str()));
}
TEST_F(PlotTaskTest, DORIRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.dori = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_DORI.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_DORI.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print DORI to file.
 */
TEST_F(PlotTaskTest, SignedDensityRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.signedDensity = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_signedDensity.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_signedDensity.cube").c_str()));
}
TEST_F(PlotTaskTest, SignedDensityRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.signedDensity = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_signedDensity.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_signedDensity.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print density to file.
 */
TEST_F(PlotTaskTest, DensityUnrestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto task = PlotTask<UNRESTRICTED>({systemController}, {});
  task.settings.density = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_alpha_Density.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_alpha_Density.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_beta_Density.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_beta_Density.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_TotalDensity.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_TotalDensity.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SpinDensity.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SpinDensity.cube").c_str()));
}
TEST_F(PlotTaskTest, DensityUnrestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto task = PlotTask<UNRESTRICTED>({systemController}, {});
  task.settings.density = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_alpha_Density.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_alpha_Density.dat").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_beta_Density.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_beta_Density.dat").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_TotalDensity.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_TotalDensity.dat").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SpinDensity.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_SpinDensity.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print MOs to file.
 */
TEST_F(PlotTaskTest, MOsRestrictedCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.allOrbitals = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.cube").c_str()));
}
TEST_F(PlotTaskTest, MOsRestrictedPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.allOrbitals = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.dat").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print MOs to file.
 */
TEST_F(PlotTaskTest, MOsRestricted_partsCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.orbitals = {1};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.cube").c_str()));
}
TEST_F(PlotTaskTest, MOsRestricted_partsPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.orbitals = {1};
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.dat").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print MOs to file.
 */
TEST_F(PlotTaskTest, MOsRestricted_nonExistentCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.orbitals = {1, 2, 3};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.cube").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.cube").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.cube").c_str()));
}
TEST_F(PlotTaskTest, MOsRestricted_nonExistentPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.orbitals = {1, 2, 3};
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO1.dat").c_str()));
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO2.dat").c_str()));
  EXPECT_FALSE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_MO3.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print ESP to file.
 */
TEST_F(PlotTaskTest, ESPCube) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.electrostaticPot = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ESP.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ESP.cube").c_str()));
}
TEST_F(PlotTaskTest, ESPPlane) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  task.settings.electrostaticPot = true;
  task.settings.p1 = {0.0, 0.0, 0.0};
  task.settings.p2 = {0.0, 1.0, 0.0};
  task.settings.p3 = {0.0, 0.0, 1.0};
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ESP.dat").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS_ESP.dat").c_str()));
}

/**
 * @test
 * @brief Tests PlotTask.h/.cpp: Print Supersystem NTOs and MOs in the restricted case, as well as transition, particle
 * and hole densities.
 */
TEST_F(PlotTaskTest, rNTOs_plusMOs) {
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  std::vector<std::shared_ptr<SystemController>> active = {systemController};
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  lrscf.settings.nEigen = 5;
  lrscf.run();
  task.settings.ntos = true;
  task.settings.excitations = {2, 3, 5};
  task.settings.orbitals = {13, 62};
  task.run();
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + systemController->getSystemName() + "_MO13.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + systemController->getSystemName() + "_MO62.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO3/7_occ.cube").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO3/10_virt.cube").c_str()));
  for (unsigned iExc : task.settings.excitations) {
    EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO" + std::to_string(iExc) + "/README").c_str()));
    EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO" + std::to_string(iExc) + "/8_occ.cube").c_str()));
    EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO" + std::to_string(iExc) + "/9_virt.cube").c_str()));
    EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS/NTO" + std::to_string(iExc)).c_str()));
    for (unsigned holeparticle = 0; holeparticle < 3; holeparticle++) {
      std::string filename =
          (holeparticle == 0) ? "_HoleDensity_" : (holeparticle == 1 ? "_ParticleDensity_" : "_TransitionDensity_");
      EXPECT_TRUE(fileExists((systemController->getSystemPath() + systemController->getSystemName() + filename +
                              std::to_string(iExc) + ".cube")
                                 .c_str()));
      EXPECT_EQ(0, std::remove((systemController->getSystemPath() + systemController->getSystemName() + filename +
                                std::to_string(iExc) + ".cube")
                                   .c_str()));
    }
  }
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "/NTOS").c_str()));
}

} /* namespace Serenity */
