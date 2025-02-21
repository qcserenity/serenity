/**
 * @file LRSCFTaskSubTDDFT_test.cpp
 *
 * @date May 31, 2018
 * @author J. Toelle
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
#include "data/ElectronicStructure.h"
#include "io/HDF5.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "tasks/VirtualOrbitalSpaceSelectionTask.h" //To be tested.
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFTask.h/.cpp .
 */
class LRSCFTaskSubTDDFTTest : public ::testing::Test {
 protected:
  LRSCFTaskSubTDDFTTest() {
  }

  virtual ~LRSCFTaskSubTDDFTTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_LLP91) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86, true);

  auto task = FreezeAndThawTask<RESTRICTED>({act, env});
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.maxCycles = 3;
  task.run();

  // Uncoupled Subsystem A.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.nEigen = 3;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscfA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  // Uncoupled Subsystem B.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.nEigen = 3;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscfB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  // Coupled Subsystems A and B.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  // Test partial response construction along ..
  lrscfAB.settings.partialResponseConstruction = true;
  lrscfAB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBA({act, env}, {});
  lrscfBA.settings.method = Options::LR_METHOD::TDA;
  lrscfBA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfBA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscfBA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfBA.settings.grid.smallGridAccuracy = 7;
  lrscfBA.settings.grid.accuracy = 7;
  // Test partial response construction along ..
  lrscfBA.settings.partialResponseConstruction = false;
  lrscfBA.run();

  // Serenity Feb 2021.
  Eigen::VectorXd excActUncoupled(3);
  excActUncoupled << 0.4946720, 0.7543262, 1.1714099;
  Eigen::VectorXd excEnvUncoupled(3);
  excEnvUncoupled << 0.7162539, 0.8659201, 1.6219693;
  Eigen::VectorXd excCoupled(6);
  excCoupled << 0.4730031, 0.6331291, 0.8006612, 0.8541067, 1.1184357, 1.7452153;

  // Compare ground-state energies.
  EXPECT_NEAR(-1.8095641208, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  EXPECT_NEAR(-1.8031312281, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  // Compare excitation energies.
  double maxDiffActUncoupled = (excActUncoupled - lrscfA.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffEnvUncoupled = (excEnvUncoupled - lrscfB.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffCoupled = (excCoupled - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffPartialConstructionVSregular =
      (lrscfBA.getTransitions().col(0) - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff();

  EXPECT_LE(maxDiffActUncoupled, 1e-6);
  EXPECT_LE(maxDiffEnvUncoupled, 1e-6);
  EXPECT_LE(maxDiffCoupled, 1e-6);
  EXPECT_LE(maxDiffPartialConstructionVSregular, 1e-8);

  std::remove((act->getSystemName() + "_" + env->getSystemName() + "_FDEcMatrix.txt").c_str());
  std::remove((act->getSystemName() + "/" + env->getSystemName() + ".TDACoupling.txt").c_str());
  std::remove((env->getSystemName() + "/" + act->getSystemName() + ".TDACoupling.txt").c_str());
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_LLP91_customFunctional) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86, true);

  auto task = FreezeAndThawTask<RESTRICTED>({act, env});
  task.settings.embedding.customNaddXCFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                               BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  task.settings.embedding.customNaddXCFunc.mixingFactors = {1.0, 1.0};
  task.settings.embedding.customNaddXCFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  task.settings.embedding.customNaddKinFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::K_LLP};
  task.settings.embedding.customNaddKinFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.maxCycles = 3;
  task.run();

  // Uncoupled Subsystem A.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.nEigen = 3;
  lrscfA.settings.embedding.customNaddXCFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                                 BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  lrscfA.settings.embedding.customNaddXCFunc.mixingFactors = {1.0, 1.0};
  lrscfA.settings.embedding.customNaddXCFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfA.settings.embedding.customNaddKinFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::K_LLP};
  lrscfA.settings.embedding.customNaddKinFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  // Uncoupled Subsystem B.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.nEigen = 3;
  lrscfB.settings.embedding.customNaddXCFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                                 BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  lrscfB.settings.embedding.customNaddXCFunc.mixingFactors = {1.0, 1.0};
  lrscfB.settings.embedding.customNaddXCFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfB.settings.embedding.customNaddKinFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::K_LLP};
  lrscfB.settings.embedding.customNaddKinFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  // Coupled Subsystems A and B.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.embedding.customNaddXCFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                                  BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  lrscfAB.settings.embedding.customNaddXCFunc.mixingFactors = {1.0, 1.0};
  lrscfAB.settings.embedding.customNaddXCFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfAB.settings.embedding.customNaddKinFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::K_LLP};
  lrscfAB.settings.embedding.customNaddKinFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  // Test partial response construction along ..
  lrscfAB.settings.partialResponseConstruction = true;
  lrscfAB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBA({act, env}, {});
  lrscfBA.settings.method = Options::LR_METHOD::TDA;
  lrscfBA.settings.embedding.customNaddXCFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                                  BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  lrscfBA.settings.embedding.customNaddXCFunc.mixingFactors = {1.0, 1.0};
  lrscfBA.settings.embedding.customNaddXCFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfBA.settings.embedding.customNaddKinFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::K_LLP};
  lrscfBA.settings.embedding.customNaddKinFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfBA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfBA.settings.grid.smallGridAccuracy = 7;
  lrscfBA.settings.grid.accuracy = 7;
  // Test partial response construction along ..
  lrscfBA.settings.partialResponseConstruction = false;
  lrscfBA.run();

  // Serenity Feb 2021.
  Eigen::VectorXd excActUncoupled(3);
  excActUncoupled << 0.4946720, 0.7543262, 1.1714099;
  Eigen::VectorXd excEnvUncoupled(3);
  excEnvUncoupled << 0.7162539, 0.8659201, 1.6219693;
  Eigen::VectorXd excCoupled(6);
  excCoupled << 0.4730031, 0.6331291, 0.8006612, 0.8541067, 1.1184357, 1.7452153;

  // Compare ground-state energies.
  EXPECT_NEAR(-1.8095641208, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  EXPECT_NEAR(-1.8031312281, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  // Compare excitation energies.
  double maxDiffActUncoupled = (excActUncoupled - lrscfA.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffEnvUncoupled = (excEnvUncoupled - lrscfB.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffCoupled = (excCoupled - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff();
  double maxDiffPartialConstructionVSregular =
      (lrscfBA.getTransitions().col(0) - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff();

  EXPECT_LE(maxDiffActUncoupled, 1e-6);
  EXPECT_LE(maxDiffEnvUncoupled, 1e-6);
  EXPECT_LE(maxDiffCoupled, 1e-6);
  EXPECT_LE(maxDiffPartialConstructionVSregular, 1e-8);

  std::remove((act->getSystemName() + "_" + env->getSystemName() + "_FDEcMatrix.txt").c_str());
  std::remove((act->getSystemName() + "/" + env->getSystemName() + ".TDACoupling.txt").c_str());
  std::remove((env->getSystemName() + "/" + act->getSystemName() + ".TDACoupling.txt").c_str());
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_subsystemGrid) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86, true);

  auto task = FreezeAndThawTask<RESTRICTED>({act, env});
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.maxCycles = 3;
  task.run();

  // Uncoupled Subsystem A.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.nEigen = 3;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscfA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.settings.subsystemgrid = {1};
  lrscfA.run();

  // Serenity Feb 2021.
  Eigen::VectorXd excActUncoupled(3);
  excActUncoupled << 0.4946720, 0.7543262, 1.1714099;
  Eigen::VectorXd excEnvUncoupled(3);
  excEnvUncoupled << 0.7270414, 0.8734773, 1.6619003;
  Eigen::VectorXd excCoupled(6);
  excCoupled << 0.4732598, 0.6351009, 0.8004093, 0.8541049, 1.1200301, 1.7431173;

  // Compare ground-state energies.
  EXPECT_NEAR(-1.8095641208, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);
  EXPECT_NEAR(-1.8031312281, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  // Compare excitation energies from "full grid" to "grid only for susbystem act".
  double maxDiffActUncoupled = (excActUncoupled - lrscfA.getTransitions().col(0)).cwiseAbs().maxCoeff();

  EXPECT_LE(maxDiffActUncoupled, 1e-6);
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Level) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 8e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Level_Triplet) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.settings.triplet = true;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.settings.triplet = true;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.triplet = true;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.3137398, 0.4916258, 0.7148767, 0.8541922, 1.0386668, 1.1007225, 1.1614874, 1.2432037,
      1.6542859, 1.7057426, 2.9531630, 3.4006584;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 8e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Huz) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();

  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Fermi_Shifted_Huz) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();

  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Hof) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();

  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();

  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.method = Options::LR_METHOD::TDA;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDDFT_Level) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP, true);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  task.run();

  EXPECT_NEAR(-1.8049132138, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  task2.run();

  EXPECT_NEAR(-1.8049132138, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.partialResponseConstruction = true;

  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.41363795, 0.58140045, 0.75028406, 0.95883377, 1.12931793, 1.19638106, 1.32308888,
      1.36862945, 1.74765800, 1.88834378, 3.11266512, 3.63866384;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 4e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDDFT_Level_Unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP, true);
  auto env =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP, true);

  act->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  env->setSCFMode(Options::SCF_MODES::UNRESTRICTED);

  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  task.run();

  EXPECT_NEAR(-1.8049132138, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<UNRESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  task2.run();

  EXPECT_NEAR(-1.8049132138, env->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 12;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 12;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 24;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.run();

  auto excitationsActUncoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(24);
  excitationEnergy_reference << 0.3158297, 0.4136380, 0.4981369, 0.5814007, 0.6960251, 0.7502841, 0.8674703, 0.9588338,
      1.0585691, 1.0970463, 1.1293179, 1.1647703, 1.1963811, 1.2487130, 1.3230889, 1.3686295, 1.6684749, 1.7343388,
      1.7476579, 1.8883437, 2.9820763, 3.1126655, 3.4178812, 3.6386643;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 4e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Huzinaga_OrbitalSelection) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  // Sys LE
  auto act_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings_act_LE = act_LE->getSettings();
  settings_act_LE.name = act_LE->getSettings().name + "LE";
  act_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86,
                                                              settings_act_LE);
  // Sys CT
  auto act_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings_act_CT = act_CT->getSettings();
  settings_act_CT.name = act_CT->getSettings().name + "CT";
  act_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86,
                                                              settings_act_CT);

  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  // Sys LE
  auto env_LE =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings settings_env_LE = env_LE->getSettings();
  settings_env_LE.name = env_LE->getSettings().name + "LE";
  env_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86,
                                                              settings_env_LE);
  // Sys CT
  auto env_CT =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings settings_env_CT = env_CT->getSettings();
  settings_env_CT.name = env_CT->getSettings().name + "CT";
  env_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86,
                                                              settings_env_CT);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossALE({act, act_LE}, {env});
  vossALE.settings.excludeProjection = true;
  vossALE.settings.localizedVirtualorbitals = true;
  vossALE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossALE.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  vossALE.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossACT({act, act_CT}, {env});
  vossACT.settings.excludeProjection = true;
  vossACT.settings.localizedEnvVirtualorbitals = true;
  vossACT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossACT.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  vossACT.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossBLE({env, env_LE}, {act});
  vossBLE.settings.excludeProjection = true;
  vossBLE.settings.localizedVirtualorbitals = true;
  vossBLE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossBLE.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  vossBLE.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossBCT({env, env_CT}, {act});
  vossBCT.settings.excludeProjection = true;
  vossBCT.settings.localizedEnvVirtualorbitals = true;
  vossBCT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossBCT.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  vossBCT.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfALE({act_LE}, {env});
  lrscfALE.settings.method = Options::LR_METHOD::TDA;
  lrscfALE.settings.nEigen = 6;
  lrscfALE.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfALE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfALE.settings.grid.smallGridAccuracy = 7;
  lrscfALE.settings.grid.accuracy = 7;
  lrscfALE.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfACT({act_CT}, {env});
  lrscfACT.settings.method = Options::LR_METHOD::TDA;
  lrscfACT.settings.nEigen = 6;
  lrscfACT.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfACT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfACT.settings.grid.smallGridAccuracy = 7;
  lrscfACT.settings.grid.accuracy = 7;
  lrscfACT.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBLE({env_LE}, {act});
  lrscfBLE.settings.method = Options::LR_METHOD::TDA;
  lrscfBLE.settings.nEigen = 6;
  lrscfBLE.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfBLE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfBLE.settings.grid.smallGridAccuracy = 7;
  lrscfBLE.settings.grid.accuracy = 7;
  lrscfBLE.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBCT({env_CT}, {act});
  lrscfBCT.settings.method = Options::LR_METHOD::TDA;
  lrscfBCT.settings.nEigen = 6;
  lrscfBCT.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfBCT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfBCT.settings.grid.smallGridAccuracy = 7;
  lrscfBCT.settings.grid.accuracy = 7;
  lrscfBCT.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfCoupling({act_LE, act_CT, env_LE, env_CT}, {});
  lrscfCoupling.settings.method = Options::LR_METHOD::TDA;
  lrscfCoupling.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfCoupling.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfCoupling.settings.grid.smallGridAccuracy = 7;
  lrscfCoupling.settings.grid.accuracy = 7;
  lrscfCoupling.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  lrscfCoupling.settings.samedensity = {1, 1, 3, 3};
  lrscfCoupling.run();

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act_LE, lrscfCoupling.settings);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsActCoupled)(iExc), 1e-6);
  }

  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act_LE);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act_CT);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env_LE);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env_CT);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, sTDA_Level_OrbitalSelection) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  // Sys LE
  auto act_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings_act_LE = act_LE->getSettings();
  settings_act_LE.name = act_LE->getSettings().name + "LE";
  act_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86,
                                                              settings_act_LE);
  // Sys CT
  auto act_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  Settings settings_act_CT = act_CT->getSettings();
  settings_act_CT.name = act_CT->getSettings().name + "CT";
  act_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86,
                                                              settings_act_CT);

  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  // Sys LE
  auto env_LE =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings settings_env_LE = env_LE->getSettings();
  settings_env_LE.name = env_LE->getSettings().name + "LE";
  env_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86,
                                                              settings_env_LE);
  // Sys CT
  auto env_CT =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings settings_env_CT = env_CT->getSettings();
  settings_env_CT.name = env_CT->getSettings().name + "CT";
  env_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86,
                                                              settings_env_CT);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task2.run();
  // From Serenity 02.08.2019
  EXPECT_NEAR(-1.8155753808, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossALE({act, act_LE}, {env});
  vossALE.settings.excludeProjection = true;
  vossALE.settings.localizedVirtualorbitals = true;
  vossALE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossALE.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossALE.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossACT({act, act_CT}, {env});
  vossACT.settings.excludeProjection = true;
  vossACT.settings.localizedEnvVirtualorbitals = true;
  vossACT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossACT.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossACT.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossBLE({env, env_LE}, {act});
  vossBLE.settings.excludeProjection = true;
  vossBLE.settings.localizedVirtualorbitals = true;
  vossBLE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossBLE.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossBLE.run();
  VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED> vossBCT({env, env_CT}, {act});
  vossBCT.settings.excludeProjection = true;
  vossBCT.settings.localizedEnvVirtualorbitals = true;
  vossBCT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossBCT.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  vossBCT.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfALE({act_LE}, {env});
  lrscfALE.settings.method = Options::LR_METHOD::TDA;
  lrscfALE.settings.nEigen = 6;
  lrscfALE.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfALE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfALE.settings.grid.smallGridAccuracy = 7;
  lrscfALE.settings.grid.accuracy = 7;
  lrscfALE.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfACT({act_CT}, {env});
  lrscfACT.settings.method = Options::LR_METHOD::TDA;
  lrscfACT.settings.nEigen = 6;
  lrscfACT.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfACT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfACT.settings.grid.smallGridAccuracy = 7;
  lrscfACT.settings.grid.accuracy = 7;
  lrscfACT.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBLE({env_LE}, {act});
  lrscfBLE.settings.method = Options::LR_METHOD::TDA;
  lrscfBLE.settings.nEigen = 6;
  lrscfBLE.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfBLE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfBLE.settings.grid.smallGridAccuracy = 7;
  lrscfBLE.settings.grid.accuracy = 7;
  lrscfBLE.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfBCT({env_CT}, {act});
  lrscfBCT.settings.method = Options::LR_METHOD::TDA;
  lrscfBCT.settings.nEigen = 6;
  lrscfBCT.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfBCT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfBCT.settings.grid.smallGridAccuracy = 7;
  lrscfBCT.settings.grid.accuracy = 7;
  lrscfBCT.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfCoupling({act_LE, act_CT, env_LE, env_CT}, {});
  lrscfCoupling.settings.method = Options::LR_METHOD::TDA;
  lrscfCoupling.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfCoupling.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfCoupling.settings.grid.smallGridAccuracy = 7;
  lrscfCoupling.settings.grid.accuracy = 7;
  lrscfCoupling.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfCoupling.settings.samedensity = {1, 1, 3, 3};
  lrscfCoupling.run();

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act_LE, lrscfCoupling.settings);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.408893, 0.574051, 0.758401, 0.944049, 1.102164, 1.179569, 1.318655, 1.372784,
      1.722485, 1.875651, 3.086483, 3.627045;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsActCoupled)(iExc), 1e-6);
  }

  std::string name = "TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act_LE);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act_CT);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env_LE);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env_CT);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, rTDDFT_Exact_Approx_Embedding) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;

  auto act1_act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterDim_DFT, true);
  Settings settings_act1_act2 = act1_act2->getSettings();
  settings_act1_act2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  act1_act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterDim_DFT, settings_act1_act2);

  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT, true);
  Settings settings_act1 = act1->getSettings();
  settings_act1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT, settings_act1);

  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT, true);
  Settings settings_act2 = act2->getSettings();
  settings_act2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT, settings_act2);

  auto env3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT, true);
  Settings settings_env3 = env3->getSettings();
  settings_env3.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  env3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT, settings_env3);

  act1_act2->setSCFMode(SPIN);
  act1->setSCFMode(SPIN);
  act2->setSCFMode(SPIN);
  env3->setSCFMode(SPIN);
  auto task = FreezeAndThawTask<SPIN>({act1_act2, env3}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  // Perform FaT with combined 1+2 and 3 calculation
  task.run();
  // LRSCF A and B
  LRSCFTask<SPIN> lrscfAB({act1_act2}, {env3});
  lrscfAB.settings.nEigen = 5;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  lrscfAB.run();
  auto lrscfContrCombined = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act1_act2, lrscfAB.settings);
  auto excitationsABCombined = lrscfContrCombined->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  // Perform Localization
  auto localization = LocalizationTask(act1_act2);
  localization.run();
  // Splitting
  auto splitting = SystemSplittingTask<SPIN>(act1_act2, {act1, act2});
  splitting.run();
  // Run additional FaT with 1 cycle to relax exactly treated systems
  auto task2 = FDETask<SPIN>({act1}, {act2, env3});
  task2.settings.maxCycles = 1;
  task2.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task2.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::PBE0,
                                             CompositeFunctionals::XCFUNCTIONALS::PBE};
  task2.run();
  auto task3 = FDETask<SPIN>({act2}, {act1, env3});
  task3.settings.maxCycles = 1;
  task3.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task3.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task3.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::PBE0,
                                             CompositeFunctionals::XCFUNCTIONALS::PBE};
  task3.run();
  VirtualOrbitalSpaceSelectionTask<SPIN> vossA({act1}, {act2});
  vossA.settings.excludeProjection = true;
  vossA.run();
  LRSCFTask<SPIN> lrscfA({act1}, {act2, env3});
  lrscfA.settings.nEigen = 5;
  lrscfA.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfA.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::PBE0,
                                              CompositeFunctionals::XCFUNCTIONALS::PBE};
  lrscfA.run();
  VirtualOrbitalSpaceSelectionTask<SPIN> vossB({act2}, {act1});
  vossB.settings.excludeProjection = true;
  vossB.run();
  LRSCFTask<SPIN> lrscfB({act2}, {act1, env3});
  lrscfB.settings.nEigen = 5;
  lrscfB.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfB.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::PBE0,
                                              CompositeFunctionals::XCFUNCTIONALS::PBE};
  lrscfB.run();

  LRSCFTask<SPIN> lrscfAB_New({act1, act2}, {env3});
  lrscfAB_New.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                      Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                      Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfAB_New.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfAB_New.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::PBE0,
                                                   CompositeFunctionals::XCFUNCTIONALS::PBE};
  lrscfAB_New.settings.fullFDEc = true;
  lrscfAB_New.run();
  auto lrscfContrSeperated1 = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act1, lrscfAB_New.settings);
  auto excitationsAB = lrscfContrSeperated1->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  for (unsigned int iExc = 0; iExc < (*excitationsABCombined).size(); iExc++) {
    EXPECT_NEAR((*excitationsAB)(iExc), (*excitationsABCombined)(iExc), 1.0e-5);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1_act2);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, uTDDFT_Exact_Approx_Embedding) {
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
  // Perform FaT with combined 1+2 and 3 calculation
  task.run();
  // LRSCF A and B
  LRSCFTask<SPIN> lrscfAB({act1_act2}, {env3});
  lrscfAB.settings.nEigen = 999;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.run();
  auto lrscfContrCombined = std::make_shared<LRSCFController<SPIN>>(act1_act2, lrscfAB.settings);
  auto excitationsABCombined = lrscfContrCombined->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  // Perform Localization
  auto localization = LocalizationTask(act1_act2);
  localization.run();
  // Splitting
  auto splitting = SystemSplittingTask<SPIN>(act1_act2, {act1, act2});
  splitting.run();
  // Run additional FaT with 1 cycle to relax exactly treated systems
  auto task2 = FreezeAndThawTask<SPIN>({act1, act2}, {env3});
  task2.settings.maxCycles = 1;
  task2.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task2.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                             CompositeFunctionals::XCFUNCTIONALS::BP86};
  task2.run();
  VirtualOrbitalSpaceSelectionTask<SPIN> vossA({act1}, {act2});
  vossA.settings.excludeProjection = true;
  vossA.run();
  LRSCFTask<SPIN> lrscfA({act1}, {act2, env3});
  lrscfA.settings.nEigen = 999;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfA.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                              CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfA.run();
  VirtualOrbitalSpaceSelectionTask<SPIN> vossB({act2}, {act1});
  vossB.settings.excludeProjection = true;
  vossB.run();
  LRSCFTask<SPIN> lrscfB({act2}, {act1, env3});
  lrscfB.settings.nEigen = 999;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfB.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                              CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfB.run();

  LRSCFTask<SPIN> lrscfAB_New({act1, act2}, {env3});
  lrscfAB_New.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                      Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                      Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfAB_New.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfAB_New.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                                   CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfAB_New.settings.fullFDEc = true;
  lrscfAB_New.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB_New.run();
  auto lrscfContrSeperated1 = std::make_shared<LRSCFController<SPIN>>(act1, lrscfAB_New.settings);
  auto excitationsAB = lrscfContrSeperated1->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  for (unsigned int iExc = 0; iExc < (*excitationsAB).size(); iExc++) {
    EXPECT_NEAR((*excitationsAB)(iExc), (*excitationsABCombined)(iExc), 3e-5);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1_act2);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, uTDDFT_Exact_Approx_Embedding_VirtualOrbLocalization) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;
  auto act1_act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_6_31Gs_BP86, true);
  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, true);
  // LE
  auto act1_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, true);
  Settings settings_act_LE = act1_LE->getSettings();
  settings_act_LE.name = act1_LE->getSettings().name + "LE";
  act1_LE = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, settings_act_LE);
  // CT
  auto act1_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, true);
  Settings settings_act_CT = act1_CT->getSettings();
  settings_act_CT.name = act1_CT->getSettings().name + "CT";
  act1_CT = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_1_6_31Gs_BP86, settings_act_CT);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_2_6_31Gs_BP86, true);
  auto env3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_3_6_31Gs_BP86, true);

  act1_act2->setSCFMode(SPIN);
  act1->setSCFMode(SPIN);
  act1_LE->setSCFMode(SPIN);
  act1_CT->setSCFMode(SPIN);
  act2->setSCFMode(SPIN);
  env3->setSCFMode(SPIN);
  auto task = FreezeAndThawTask<SPIN>({act1_act2, env3}, {});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  // Perform FaT with combined 1+2 and 3 calculation
  task.run();

  // Perform Localization
  auto localization = LocalizationTask(act1_act2);
  localization.run();
  // Splitting
  auto splitting = SystemSplittingTask<SPIN>(act1_act2, {act1, act2});
  splitting.run();
  // Run additional FaT with 1 cycle to relax exactly treated systems
  auto task2 = FreezeAndThawTask<SPIN>({act1, act2}, {env3});
  task2.settings.maxCycles = 1;
  task2.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task2.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                             CompositeFunctionals::XCFUNCTIONALS::BP86};
  task2.run();

  VirtualOrbitalSpaceSelectionTask<SPIN> vossA({act1}, {act2});
  vossA.settings.excludeProjection = true;
  vossA.run();

  LRSCFTask<SPIN> lrscfA({act1}, {act2, env3});
  lrscfA.settings.nEigen = 999;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                 Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfA.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                              CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfA.run();

  VirtualOrbitalSpaceSelectionTask<SPIN> vossALE({act1, act1_LE}, {act2, env3});
  vossALE.settings.localizedVirtualorbitals = true;
  vossALE.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossALE.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                  Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                  Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  vossALE.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  vossALE.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                               CompositeFunctionals::XCFUNCTIONALS::BP86};
  vossALE.run();

  VirtualOrbitalSpaceSelectionTask<SPIN> vossACT({act1, act1_CT}, {act2, env3});
  vossACT.settings.localizedEnvVirtualorbitals = true;
  vossACT.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  vossACT.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                  Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                  Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  vossACT.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  vossACT.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                               CompositeFunctionals::XCFUNCTIONALS::BP86};
  vossACT.run();

  LRSCFTask<SPIN> lrscfALE({act1_LE}, {act2, env3});
  lrscfALE.settings.nEigen = 999;
  lrscfALE.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfALE.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                   Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                   Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfALE.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfALE.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                                CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfALE.run();

  LRSCFTask<SPIN> lrscfACT({act1_CT}, {act2, env3});
  lrscfACT.settings.nEigen = 999;
  lrscfACT.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfACT.settings.embedding.embeddingModeList = {Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                   Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
                                                   Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfACT.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfACT.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                                CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfACT.run();

  LRSCFTask<SPIN> lrscfCoupling({act1_CT, act1_LE}, {act2, env3});
  lrscfCoupling.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfCoupling.settings.embedding.embeddingModeList = {
      Options::KIN_EMBEDDING_MODES::LEVELSHIFT, Options::KIN_EMBEDDING_MODES::LEVELSHIFT,
      Options::KIN_EMBEDDING_MODES::LEVELSHIFT, Options::KIN_EMBEDDING_MODES::NADD_FUNC};
  lrscfCoupling.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  lrscfCoupling.settings.embedding.naddXCFuncList = {CompositeFunctionals::XCFUNCTIONALS::BP86,
                                                     CompositeFunctionals::XCFUNCTIONALS::BP86};
  lrscfCoupling.settings.samedensity = {1, 1, 3, 4};
  lrscfCoupling.run();

  auto lrscfContrACoupled = std::make_shared<LRSCFController<SPIN>>(act1, lrscfCoupling.settings);
  auto excitationsCouppledLocalized_A = lrscfContrACoupled->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);

  auto lrscfContrAUncoupled = std::make_shared<LRSCFController<SPIN>>(act1, lrscfA.settings);
  auto excitationsUncoupled_A = lrscfContrAUncoupled->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsCouppledLocalized_A).size(); iExc++) {
    EXPECT_NEAR((*excitationsUncoupled_A)(iExc), (*excitationsCouppledLocalized_A)(iExc), 1e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1_act2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1_LE);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1_CT);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env3);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskSubTDDFTTest, MRICoulombResponse) {
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon1_3_DFT, true);
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon2_3_DFT, true);
  auto sys3 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMon3_3_DFT, true);

  auto task = FreezeAndThawTask<RESTRICTED>({sys1, sys2, sys3});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  task.settings.maxCycles = 1;
  task.run();

  LRSCFTask<RESTRICTED> lrscf1({sys1, sys2, sys3});
  lrscf1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscf1.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscf1.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscf1.settings.nEigen = 0;
  lrscf1.settings.frequencies = {2};

  LRSCFTask<RESTRICTED> lrscf2({sys1, sys2, sys3});
  lrscf2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscf2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscf2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  lrscf2.settings.nEigen = 0;
  lrscf2.settings.frequencies = {2};
  lrscf2.settings.approxCoulomb = {2, 3};

  // Vanilla run.
  lrscf1.run();
  auto results1 = lrscf1.getProperties();
  double pol1 = 1.0 / 3.0 * std::get<1>(lrscf1.getProperties()[0]).trace();
  double ord1 = -1.0 / 3.0 / (lrscf1.settings.frequencies[0]) * std::get<2>(lrscf1.getProperties()[0]).trace();

  // Approximate run: one interaction exactly, one interaction via MRI, and one interaction via Grimme integrals.
  lrscf2.run();
  auto results2 = lrscf2.getProperties();
  double pol2 = 1.0 / 3.0 * std::get<1>(lrscf2.getProperties()[0]).trace();
  double ord2 = -1.0 / 3.0 / (lrscf2.settings.frequencies[0]) * std::get<2>(lrscf2.getProperties()[0]).trace();

  // Compare approximate vs "exact".
  EXPECT_NEAR(pol1, pol2, 2E-2);
  EXPECT_NEAR(ord1, ord2, 4E-4);

  // Compare approximate vs old value (April 2023).
  EXPECT_NEAR(pol2, 16.3785521496, 1E-6);
  EXPECT_NEAR(ord2, 0.0616015787, 1E-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys3);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
