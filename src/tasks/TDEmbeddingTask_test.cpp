/**
 * @file TDEmbeddingTask_test.cpp
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
#include "tasks/TDEmbeddingTask.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/CoupledClusterTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemAdditionTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
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
 *        Note: This tests print-level adjustments for the task, too.
 */
TEST_F(TDEmbeddingTaskTest, restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.parseGeneralSettings();
  task.run();
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  task.parseGeneralSettings();
  EXPECT_NEAR(-1.8155882075161134, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        with different basis sets in environment and active system.
 *        Note: This tests print-level adjustments for the task, too.
 */
TEST_F(TDEmbeddingTaskTest, restrictedDoubleHybridBasisSetChange) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  task.parseGeneralSettings();
  task.run();
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  task.parseGeneralSettings();
  task.settings.mp2Type = Options::MP2_TYPES::DF;
  // Changed (12.02.2020) due to separate evaluation of Coulomb and XX
  EXPECT_NEAR(-1.78715623152954, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        with different basis sets in environment and active system.
 *        Note: This tests print-level adjustments for the task, too.
 */
TEST_F(TDEmbeddingTaskTest, restrictedDoubleHybrid_NORI) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);

  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  Settings settingsEnv = env->getSettings();
  settingsEnv.basis.densFitJ = Options::DENS_FITS::NONE;
  settingsEnv.basis.densFitK = Options::DENS_FITS::NONE;
  settingsEnv.basis.densFitLRK = Options::DENS_FITS::NONE;
  settingsEnv.basis.densFitCorr = Options::DENS_FITS::NONE;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settingsEnv);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.parseGeneralSettings();
  task.settings.mp2Type = Options::MP2_TYPES::AO;
  task.run();
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  task.parseGeneralSettings();
  // The difference between RI and no-fitting is around 1.8e-5.
  EXPECT_NEAR(-1.7871348696588256, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        with different basis sets in environment and active system.
 *        Note: This tests print-level adjustments for the task, too.
 */
TEST_F(TDEmbeddingTaskTest, restrictedDoubleHybrid_LOCAL) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.parseGeneralSettings();
  task.run();
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  task.parseGeneralSettings();
  // The difference between RI and no-fitting is around 1.8e-5.
  EXPECT_NEAR(-1.7871413196196728, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);

  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        with different basis sets in environment and active system.
 *        Note: This tests print-level adjustments for the task, too.
 */
TEST_F(TDEmbeddingTaskTest, restrictedCCSDinDoubleHybrid_LOCAL) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  Settings settings = act->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  Settings envSettings = env->getSettings();
  envSettings.basis.label = "DEF2-SVP";
  envSettings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, envSettings);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.parseGeneralSettings();
  task.run();
  task.generalSettings.printLevel = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  task.parseGeneralSettings();
  // The difference between RI and no-fitting is around 1.8e-5.
  double naddXCEnergy =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::FDE_NAD_XC);
  EXPECT_NEAR(-1.7095757608, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);

  CoupledClusterTask ccTask(act, {env});
  ccTask.settings.level = Options::CC_LEVEL::DLPNO_CCSD;
  ccTask.settings.lcSettings.embeddingSettings = task.settings.embedding;
  ccTask.run();
  double naddXCEnergy2 =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::FDE_NAD_XC);

  std::cout << "NaddXCEnergy " << naddXCEnergy << "  " << naddXCEnergy2 << std::endl;

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(
      env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/",
      "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy.
 */
TEST_F(TDEmbeddingTaskTest, restricted_read_supersystem) {
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/";
  }
  else {
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set." << std::endl;
    assert(false);
  }
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.load = pathToTestsResources;
  task.settings.name = "TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE";
  task.run();
  EXPECT_NEAR(-1.8155882075161134, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory("TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE/",
                                                        "TestSystem_H2_6_31Gs_SUPERSYSTEM_FDE");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::VERBOSE;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  Settings settings = env->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.grid.smallGridAccuracy = 4;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, settings);
  settings.name = act->getSystemName();
  settings.path = act->getSystemPath();
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);

  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.splitValenceAndCore = true;
  task.run();
  EXPECT_NEAR(-1.8155886122, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy with Hoffmann's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_hoffmann) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.run();
  EXPECT_NEAR(-1.8155881852554803, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy with Huzinaga's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_huzinaga) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-1.8155881852542037, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga's projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  task.run();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  EXPECT_NEAR(-1.8155881853192255, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with a fermi-shifted Huzinaga projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_fermi_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  EXPECT_NEAR(-306.6721662019, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  EXPECT_NEAR(task.settings.embedding.fermiShift, 3.430175862028e-02, 1e-7);
  std::string name = "TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with a shifted Huzinaga projection operator.
 */
TEST_F(TDEmbeddingTaskTest, restricted_shifted_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_PHENYL_DEF2_SVP_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PHENOLATE_O_DEF2_SVP_BP86, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.useFermiLevel = false;
  task.run();
  EXPECT_NEAR(-306.6721662019, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  EXPECT_NEAR(task.settings.embedding.fermiShift, 1.0, 1e-7);
  std::string name = "TestSystem_Phenolate_Phenyl_Def2_SVP_BP86+TestSystem_Phenolate_O_Def2_SVP_BP86";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga's projection operator for hybrid
 * functionals.
 */
TEST_F(TDEmbeddingTaskTest, restricted_huzinaga_hybrid) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_HYBRID, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_HYBRID, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;

  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;

  task.run();
  auto supersystem = *act + *env;
  ScfTask<Options::SCF_MODES::RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.run();
  auto superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(superEnergy, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with the level-shift operator for double hybrid
 * functionals.
 */
TEST_F(TDEmbeddingTaskTest, restricted_levelshift_doubleHybrid) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.netThreshold = 1e-4;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;

  task.settings.lcSettings.pnoThreshold = 1e-10;
  task.settings.lcSettings.mullikenThreshold = 1e-5;
  task.settings.lcSettings.orbitalToShellThreshold = 1e-5;
  task.settings.lcSettings.doiPAOThreshold = 1e-5;
  task.settings.lcSettings.useProjectedOccupiedOrbitals = true;
  task.settings.splitValenceAndCore = true;
  task.settings.mp2Type = Options::MP2_TYPES::LOCAL;

  task.run();
  // Working implementation as reference. (I am not aware that any other code has this feature or any easy way
  // to compare to a supersystem calculation).
  EXPECT_NEAR(-152.57223847898166, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 5e-5);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_hardtruncated_netpopulation_huzinaga) {
  SystemController__TEST_SUPPLY::cleanUp();
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  // HF environment
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, true);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.netThreshold = 1e-2;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-152.42678230626774, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-5);
  std::string name = "TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation and
 * subsequent restart.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_restart) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  // HF environment
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.netThreshold = 1e-4;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load = ".";
  loadSettings.path = "new";

  auto newSys = std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(), newSys->getBasisController()->getNBasisFunctions());
  std::string name = "TestSystem_WaterMonOne_6_31Gs_DFT+TestSystem_WaterMonTwo_6_31Gs";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(newSys);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_restart_basisChange) {
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS, true);
  // HF environment
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.netThreshold = 1e-3;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load = ".";
  loadSettings.path = "new";

  auto newSys = std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(), newSys->getBasisController()->getNBasisFunctions());
  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(), (unsigned int)18);

  std::string name = "TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(newSys);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_levelshift_restart_basisChange) {
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS, true);
  // HF environment
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();

  auto loadSettings(act->getSettings());
  loadSettings.load = ".";
  loadSettings.path = "new";

  auto newSys = std::make_shared<SystemController>(loadSettings);

  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(), newSys->getBasisController()->getNBasisFunctions());
  EXPECT_EQ(act->getBasisController()->getNBasisFunctions(), (unsigned int)36);
  std::string name = "TestSystem_WaterMonTwo_6_31Gs_DFT+TestSystem_WATER_DISTORTED_MINBAS";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(newSys);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + prim. net population.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_primNetPop_levelshift) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.truncationFactor = 0.5;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.run();
  EXPECT_NEAR(-1.8133458052953828, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with Huzinaga + net population truncation.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.netThreshold = 8e-2;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.run();
  EXPECT_NEAR(-1.8189900344276131, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population + non ortho. distant +
 * kin. energy functional.
 */
TEST_F(TDEmbeddingTaskTest, restricted_truncated_netpopulation_levelshift_nonortho_distant_kinetic) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.netThreshold = 1e-3;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;

  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.longRangeNaddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.borderAtomThreshold = 0.8;
  task.settings.embedding.basisFunctionRatio = 0.4;

  task.run();
  EXPECT_NEAR(-152.81569649224861, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test unrestricted energy for the Carter potential reconstruction.
 */
TEST_F(TDEmbeddingTaskTest, unrestricted_Carter) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::RECONSTRUCTION;
  task.settings.embedding.carterCycles = 5000;
  task.run();
  EXPECT_NEAR(-1.8155881842923347, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy with levelshift + net population + non ortho. distant +
 * kin. energy functional.
 */
TEST_F(TDEmbeddingTaskTest, restricted_WYLB) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::RECONSTRUCTION;
  task.settings.embedding.lbCycles = 5000;
  task.settings.noSupRec = false;
  task.run();
  EXPECT_NEAR(-1.8062078106336934, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(TDEmbeddingTaskTest, embeddingMode_NONE) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  task.run();
  EXPECT_NEAR(-2.4352547558648183, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(TDEmbeddingTaskTest, embeddingMode_NADDFUNC) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.run();
  EXPECT_NEAR(-1.7947783381675915, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-7);
  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests the TDEmbedding with ECPs.
 */
TEST_F(TDEmbeddingTaskTest, embeddingMode_Huzinaga_ECPs) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  // There is (currently) no MINAO basis for I.
  task.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();
  // Perform supersystem calculations
  auto supersystem = *act + *env;
  ScfTask<Options::SCF_MODES::RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.run();
  auto superEnergy = supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  EXPECT_NEAR(superEnergy, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-4);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(
      env->getSystemPath() + "TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/I_FREE/", "I_FREE");
  std::string name = "TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the TDEmbedding with ECPs.
 */
TEST_F(TDEmbeddingTaskTest, embeddingMode_HuzinagaTruncated_ECPs) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  // There is (currently) no MINAO basis for I.
  task.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.netThreshold = 0.0;
  task.run();
  EXPECT_NEAR(-893.6684430205, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-4);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(
      env->getSystemPath() + "TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE/I_FREE/", "I_FREE");
  std::string name = "TestSystem_HI_Def2_SVP_PBE+TestSystem_I2_Def2_SVP_PBE";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        in active and environment system.
 */
TEST_F(TDEmbeddingTaskTest, restrictedDoubleHybridinDoubleHybrid_Levelshift) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  settings = env->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.fullMP2Coupling = true;
  task.run();

  Settings supersystemSettings = env->getSettings();
  supersystemSettings.name = act->getSystemName() + "+" + env->getSystemName();
  supersystemSettings.charge = 0;
  supersystemSettings.spin = 0;
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
  SystemAdditionTask<RESTRICTED> additionTask(supersystem, {act, env});
  additionTask.settings.addOccupiedOrbitals = true;
  additionTask.run();

  ScfTask<RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  supersystemSCF.run();

  EXPECT_NEAR(supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),
              act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 5e-5);

  std::string name = "TestSystem_WaterMonOne_Def2_SVP_B2PLYP+TestSystem_WaterMonTwo_Def2_SVP";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + "TMP_Active/", "TMP_Active");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for double hybrid embedding
 *        in active and environment system.
 */
TEST_F(TDEmbeddingTaskTest, restrictedDoubleHybridinDoubleHybrid_Huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  settings = env->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.fullMP2Coupling = true;
  task.settings.useFermiLevel = false;
  task.run();

  Settings supersystemSettings = env->getSettings();
  supersystemSettings.name = act->getSystemName() + "+" + env->getSystemName();
  supersystemSettings.charge = 0;
  supersystemSettings.spin = 0;
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
  SystemAdditionTask<RESTRICTED> additionTask(supersystem, {act, env});
  additionTask.settings.addOccupiedOrbitals = true;
  additionTask.run();

  ScfTask<RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  supersystemSCF.run();

  EXPECT_NEAR(supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),
              act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 5e-5);

  std::string name = "TestSystem_WaterMonOne_Def2_SVP_B2PLYP+TestSystem_WaterMonTwo_Def2_SVP";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + "TMP_Active/", "TMP_Active");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests TDEmbeddingTask.h/.cpp: Test restricted energy for gga embedding with dispersion
 *        in active and environment system.
 */
TEST_F(TDEmbeddingTaskTest, restrictedGGADispersion) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  Settings settings = act->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.smallGridAccuracy = 4;
  settings.grid.accuracy = 4;
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  settings = env->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.smallGridAccuracy = 4;
  settings.grid.accuracy = 4;
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings);
  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJABC;
  task.settings.useFermiLevel = false;
  task.run();

  Settings supersystemSettings = env->getSettings();
  supersystemSettings.name = act->getSystemName() + "+" + env->getSystemName();
  supersystemSettings.charge = 0;
  supersystemSettings.spin = 0;
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
  SystemAdditionTask<RESTRICTED> additionTask(supersystem, {act, env});
  additionTask.settings.addOccupiedOrbitals = true;
  additionTask.run();

  ScfTask<RESTRICTED> supersystemSCF(supersystem);
  supersystemSCF.run();

  EXPECT_NEAR(supersystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),
              act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-8);

  std::string name = "TestSystem_WaterMonOne_Def2_SVP+TestSystem_WaterMonTwo_Def2_SVP";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + "TMP_Active/", "TMP_Active");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
