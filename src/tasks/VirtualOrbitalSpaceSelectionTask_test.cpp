/**
 * @file VirtualOrbitalSpaceSelectionTask_test.cpp
 *
 * @date 30.03.2020
 * @author Johannes Tölle
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
#include "tasks/VirtualOrbitalSpaceSelectionTask.h" //To be tested.
#include "basis/BasisController.h"                  //GetNBasisFunctions.
#include "data/ElectronicStructure.h"               //GetDensityMatrix.
#include "data/OrbitalController.h"                 //GetCoefficients.
#include "settings/Settings.h"
#include "system/SystemController.h" //GetElectronicStructure.
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/ScfTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class VirtualOrbitalSpaceSelectionTaskTest : public ::testing::Test {
 protected:
  VirtualOrbitalSpaceSelectionTaskTest() {
  }
  virtual ~VirtualOrbitalSpaceSelectionTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

// Most of the functionalities of this task are tested in the LRSCFTask_Test.h
TEST_F(VirtualOrbitalSpaceSelectionTaskTest, rExcludeProjection) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto nVirt = act->getNVirtualOrbitals<MODE>();

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act}, {env});
  vossA.settings.excludeProjection = true;
  vossA.run();
  auto orbEner = act->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOccEnv = env->getNOccupiedOrbitals<MODE>();
  auto nOccVirt = act->getNVirtualOrbitalsTruncated<MODE>();

  for_spin(nOccVirt, nOccEnv, nVirt) {
    EXPECT_EQ(nOccVirt_spin, nVirt_spin - nOccEnv_spin);
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(VirtualOrbitalSpaceSelectionTaskTest, uExcludeProjection) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::UNRESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<UNRESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto nVirt = act->getNVirtualOrbitals<MODE>();

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act}, {env});
  vossA.settings.excludeProjection = true;
  vossA.run();
  auto orbEner = act->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOccEnv = env->getNOccupiedOrbitals<MODE>();
  auto nOccVirt = act->getNVirtualOrbitalsTruncated<MODE>();

  for_spin(nOccVirt, nOccEnv, nVirt) {
    EXPECT_EQ(nOccVirt_spin, nVirt_spin - nOccEnv_spin);
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(VirtualOrbitalSpaceSelectionTaskTest, rVirtualOrbitalLocalization) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings = act2->getSettings();
  act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<MODE>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<MODE>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act, act2}, {env});
  vossA.settings.excludeProjection = true;
  vossA.settings.localizedVirtualorbitals = true;
  vossA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossA.run();

  auto orbEner_A = act2->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOcc = act2->getNOccupiedOrbitals<MODE>();
  auto nVirt = act2->getNVirtualOrbitalsTruncated<MODE>();

  Eigen::VectorXd ref(4);
  ref << -0.2769432308, +0.3299831617, +0.7494383145, +1.3249932782;

  for_spin(orbEner_A, nOcc, nVirt) {
    for (unsigned int counter = 0; counter < (nOcc_spin + nVirt_spin); counter++) {
      EXPECT_NEAR(ref(counter), orbEner_A_spin(counter), 1e-6);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(VirtualOrbitalSpaceSelectionTaskTest, uVirtualOrbitalLocalization) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::UNRESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings = act2->getSettings();
  act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, settings);

  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<MODE>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<MODE>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act, act2}, {env});
  vossA.settings.excludeProjection = true;
  vossA.settings.localizedVirtualorbitals = true;
  vossA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossA.run();

  auto orbEner_A = act2->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOcc = act2->getNOccupiedOrbitals<MODE>();
  auto nVirt = act2->getNVirtualOrbitalsTruncated<MODE>();
  // The same values can be used for alpha and beta spin due to degeneracy
  Eigen::VectorXd ref(4);
  ref << -0.2769432308, +0.3299831617, +0.7494383145, +1.3249932782;

  for_spin(orbEner_A, nOcc, nVirt) {
    for (unsigned int counter = 0; counter < (nOcc_spin + nVirt_spin); counter++) {
      EXPECT_NEAR(ref(counter), orbEner_A_spin(counter), 1e-6);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
//
TEST_F(VirtualOrbitalSpaceSelectionTaskTest, rCanonSelection) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings = act2->getSettings();
  act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, settings);

  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<MODE>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<MODE>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act, act2}, {env});
  vossA.settings.excludeProjection = true;
  vossA.settings.localCanonicalVirtuals = 0.51;
  vossA.run();

  auto orbEner_A = act2->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOcc = act2->getNOccupiedOrbitals<MODE>();
  auto nVirt = act2->getNVirtualOrbitalsTruncated<MODE>();

  Eigen::VectorXd ref(4);
  ref << -0.2769432308, +0.3018064565, +0.9299062933, +1.0039878814;

  for_spin(orbEner_A, nOcc, nVirt) {
    for (unsigned int counter = 0; counter < (nOcc_spin + nVirt_spin); counter++) {
      EXPECT_NEAR(ref(counter), orbEner_A_spin(counter), 1e-6);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}
//
TEST_F(VirtualOrbitalSpaceSelectionTaskTest, uCanonSelection) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::UNRESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings = act2->getSettings();
  act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, settings);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<MODE>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  auto task2 = FDETask<MODE>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<MODE>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act, act2}, {env});
  vossA.settings.excludeProjection = true;
  vossA.settings.localCanonicalVirtuals = 0.51;
  vossA.run();

  auto orbEner_A = act2->getActiveOrbitalController<MODE>()->getEigenvalues();
  auto nOcc = act2->getNOccupiedOrbitals<MODE>();
  auto nVirt = act2->getNVirtualOrbitalsTruncated<MODE>();

  Eigen::VectorXd ref(4);
  ref << -0.2769432308, +0.3018064565, +0.9299062933, +1.0039878814;

  for_spin(orbEner_A, nOcc, nVirt) {
    for (unsigned int counter = 0; counter < (nOcc_spin + nVirt_spin); counter++) {
      EXPECT_NEAR(ref(counter), orbEner_A_spin(counter), 1e-6);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(VirtualOrbitalSpaceSelectionTaskTest, rMixingOccAndVirtuals) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  auto act_env = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular, true);

  auto task = FreezeAndThawTask<MODE>({act, env}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act_env}, {act, env});
  vossA.settings.mixingOccAndVirtualOrbitals = true;
  vossA.settings.relaxation = false;
  vossA.run();

  // Reference orbital values (occ from A and virt from B)
  Eigen::VectorXd ref(4);
  ref << -0.3381540332, 0.1844602302, 0.4211059725, 1.7274564908;

  auto orbEner_A_B = act_env->getActiveOrbitalController<MODE>()->getEigenvalues();
  for_spin(orbEner_A_B) {
    for (unsigned int i = 0; i < ref.size(); i++) {
      EXPECT_NEAR(ref(i), orbEner_A_B_spin(i), 1e-6);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(VirtualOrbitalSpaceSelectionTaskTest, rMixingOccAndVirtualsRelaxation) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto MODE = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  auto act_env = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86_Supermolecular, true);

  auto task = FreezeAndThawTask<MODE>({act, env}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  VirtualOrbitalSpaceSelectionTask<MODE> vossA({act_env}, {act, env});
  vossA.settings.mixingOccAndVirtualOrbitals = true;
  vossA.settings.relaxation = true;
  vossA.run();

  // Reference orbital values (occ from A and virt from B)
  Eigen::VectorXd ref(4);
  ref << -0.442447, 0.193363, 0.340337, 1.60898;

  auto orbEner_A_B = act_env->getActiveOrbitalController<MODE>()->getEigenvalues();
  for_spin(orbEner_A_B) {
    for (unsigned int i = 0; i < ref.size(); i++) {
      EXPECT_NEAR(ref(i), orbEner_A_B_spin(i), 1e-5);
    }
  };
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86_SupersystemBasis/"
              "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86_SupersystemBasis.settings");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86_SupersystemBasis/"
              "TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86_SupersystemBasis.xyz");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86_SupersystemBasis");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/