/**
 * @file LRSCFTaskRI_test.cpp
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
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFRITask.h/.cpp .
 */
class LRSCFTaskRITest : public ::testing::Test {
 protected:
  LRSCFTaskRITest() {
  }

  virtual ~LRSCFTaskRITest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  double accuracyRICachevsNOCache = 1e-4;
  double accuracyALPHAvsBETA = 1e-6;
};

TEST_F(LRSCFTaskRITest, restricted) {
  // system
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  std::vector<std::shared_ptr<SystemController>> sysv;
  sysv.push_back(sys);

  // tddft no cache
  LRSCFTask<RESTRICTED> nocache_tddft(sysv);
  nocache_tddft.run();

  // tddft densFitK
  LRSCFTask<RESTRICTED> cache_tddft(sysv);
  cache_tddft.settings.densFitK = Options::DENS_FITS::RI;
  cache_tddft.settings.densFitLRK = Options::DENS_FITS::RI;
  cache_tddft.run();

  // compare results
  double maxDiff_TDDFT = (nocache_tddft.getTransitions() - cache_tddft.getTransitions()).cwiseAbs().maxCoeff();
  EXPECT_LE(maxDiff_TDDFT, accuracyRICachevsNOCache);

  // tda no cache
  LRSCFTask<RESTRICTED> nocache_tda(sysv);
  nocache_tda.settings.method = Options::LR_METHOD::TDA;
  nocache_tda.run();

  // tda cache
  LRSCFTask<RESTRICTED> ricache_tda(sysv);
  ricache_tda.settings.method = Options::LR_METHOD::TDA;
  ricache_tda.settings.densFitK = Options::DENS_FITS::RI;
  ricache_tda.settings.densFitLRK = Options::DENS_FITS::RI;
  ricache_tda.run();

  // compare results
  double maxDiff_TDA = (nocache_tda.getTransitions() - ricache_tda.getTransitions()).cwiseAbs().maxCoeff();
  EXPECT_LE(maxDiff_TDA, accuracyRICachevsNOCache);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskRITest, unrestricted) {
  // system
  auto sys_a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  Settings settings = sys_a->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  sys_a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  auto sys_b =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  sys_a->setSpin(1);
  sys_a->setCharge(1);
  sys_a->setSCFMode(Options::SCF_MODES::UNRESTRICTED);

  sys_b->setSpin(-1);
  sys_b->setCharge(1);
  sys_b->setSCFMode(Options::SCF_MODES::UNRESTRICTED);

  ScfTask<UNRESTRICTED> scf_a(sys_a);
  scf_a.run();

  ScfTask<UNRESTRICTED> scf_b(sys_b);
  scf_b.run();

  std::vector<std::shared_ptr<SystemController>> sysv_a;
  std::vector<std::shared_ptr<SystemController>> sysv_b;
  sysv_a.push_back(sys_a);
  sysv_b.push_back(sys_b);

  // matrices for results
  Eigen::MatrixXd ricache_tddft_a, ricache_tda_a, ricache_tddft_b, ricache_tda_b;
  Eigen::MatrixXd nocache_tddft_a, nocache_tda_a, nocache_tddft_b, nocache_tda_b;

  // separate scope
  {
    // tddft no cache
    LRSCFTask<UNRESTRICTED> nocache_tddft(sysv_a);
    nocache_tddft.run();
    nocache_tddft_a = nocache_tddft.getTransitions();

    // tddft cache
    LRSCFTask<UNRESTRICTED> ricache_tddft(sysv_a);
    ricache_tddft.settings.densFitK = Options::DENS_FITS::RI;
    ricache_tddft.settings.densFitLRK = Options::DENS_FITS::RI;
    ricache_tddft.run();
    ricache_tddft_a = ricache_tddft.getTransitions();

    // compare results
    double maxDiff_TDDFT = (nocache_tddft_a - ricache_tddft_a).cwiseAbs().maxCoeff();
    EXPECT_LE(maxDiff_TDDFT, accuracyRICachevsNOCache);

    // tda no cache
    LRSCFTask<UNRESTRICTED> nocache_tda(sysv_a);
    nocache_tda.settings.method = Options::LR_METHOD::TDA;
    nocache_tda.run();
    nocache_tda_a = nocache_tda.getTransitions();

    // tda cache
    LRSCFTask<UNRESTRICTED> ricache_tda(sysv_a);
    ricache_tda.settings.method = Options::LR_METHOD::TDA;
    ricache_tda.settings.densFitK = Options::DENS_FITS::RI;
    ricache_tda.settings.densFitLRK = Options::DENS_FITS::RI;
    ricache_tda.run();
    ricache_tda_a = ricache_tda.getTransitions();

    double maxDiff_TDA = (nocache_tda_a - ricache_tda_a).cwiseAbs().maxCoeff();
    EXPECT_LE(maxDiff_TDA, accuracyRICachevsNOCache);
  }

  // separate scope
  {
    // tddft no cache
    LRSCFTask<UNRESTRICTED> nocache_tddft(sysv_b);
    nocache_tddft.run();
    nocache_tddft_b = nocache_tddft.getTransitions();

    // tddft cache
    LRSCFTask<UNRESTRICTED> ricache_tddft(sysv_b);
    ricache_tddft.settings.densFitK = Options::DENS_FITS::RI;
    ricache_tddft.settings.densFitLRK = Options::DENS_FITS::RI;
    ricache_tddft.run();
    ricache_tddft_b = ricache_tddft.getTransitions();

    // compare results
    double maxDiff_TDDFT = (nocache_tddft_b - ricache_tddft_b).cwiseAbs().maxCoeff();
    EXPECT_LE(maxDiff_TDDFT, accuracyRICachevsNOCache);

    // tda no cache
    LRSCFTask<UNRESTRICTED> nocache_tda(sysv_b);
    nocache_tda.settings.method = Options::LR_METHOD::TDA;
    nocache_tda.run();
    nocache_tda_b = nocache_tda.getTransitions();

    // tda cache
    LRSCFTask<UNRESTRICTED> ricache_tda(sysv_b);
    ricache_tda.settings.method = Options::LR_METHOD::TDA;
    ricache_tda.settings.densFitK = Options::DENS_FITS::RI;
    ricache_tda.settings.densFitLRK = Options::DENS_FITS::RI;
    ricache_tda.run();
    ricache_tda_b = ricache_tda.getTransitions();

    double maxDiff_TDA = (nocache_tda_b - ricache_tda_b).cwiseAbs().maxCoeff();
    EXPECT_LE(maxDiff_TDA, accuracyRICachevsNOCache);
  }

  // compare results with alpha and beta excess
  double nocache_TDDFT = (nocache_tddft_a - nocache_tddft_b).cwiseAbs().maxCoeff();
  double ricache_TDDFT = (ricache_tddft_a - ricache_tddft_b).cwiseAbs().maxCoeff();
  double nocache_TDA = (nocache_tda_a - nocache_tda_b).cwiseAbs().maxCoeff();
  double ricache_TDA = (ricache_tda_a - ricache_tda_b).cwiseAbs().maxCoeff();

  EXPECT_LE(nocache_TDDFT, accuracyALPHAvsBETA);
  EXPECT_LE(ricache_TDDFT, accuracyALPHAvsBETA);
  EXPECT_LE(nocache_TDA, accuracyALPHAvsBETA);
  EXPECT_LE(ricache_TDA, accuracyALPHAvsBETA);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys_a);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys_b);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskRITest, sTDDFT_Level_densFitK_Restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings sA = act->getSettings();
  Settings sB = env->getSettings();
  sA.basis.label = "6-31G";
  sB.basis.label = "6-31G";
  sA.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  sB.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, sA);
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, sB);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  task.run();
  // From Serenity Oct 21
  EXPECT_NEAR(-1.7976759266, act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  task2.run();
  // From Serenity Oct 21
  EXPECT_NEAR(-1.7976759266, env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfA.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfA.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfA.settings.nEigen = 999;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfB.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfB.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfB.settings.nEigen = 999;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfAB.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfAB.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfAB.settings.partialResponseConstruction = true;
  lrscfAB.run();

  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference << 0.4188248, 0.5840582, 0.7474868, 0.9640923, 1.1481851, 1.2010459, 1.3219816, 1.3748426,
      1.7597483, 1.8995670, 3.1084058, 3.6524039;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskRITest, sTDDFT_Level_densFitK_Unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
  Settings sA = act->getSettings();
  Settings sB = env->getSettings();
  sA.basis.label = "6-31G";
  sB.basis.label = "6-31G";
  sA.scfMode = Options::SCF_MODES::UNRESTRICTED;
  sB.scfMode = Options::SCF_MODES::UNRESTRICTED;
  sA.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  sB.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, sA);
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, sB);

  auto task = TDEmbeddingTask<UNRESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  task.run();
  // From Serenity Oct 21
  EXPECT_NEAR(-1.7976759266, act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);

  auto task2 = FDETask<UNRESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  task2.run();
  // From Serenity Oct 21
  EXPECT_NEAR(-1.7976759266, env->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(), 1e-6);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfA.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfA.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfA.settings.nEigen = 999;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfB.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfB.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfB.settings.nEigen = 999;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfAB({act, env}, {});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(env, lrscfAB.settings);

  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::ACD;
  lrscfAB.settings.densFitK = Options::DENS_FITS::ACD;
  lrscfAB.settings.densFitLRK = Options::DENS_FITS::ACD;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
  lrscfAB.run();

  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  // From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(24);
  excitationEnergy_reference << 0.3185942, 0.4188249, 0.5065803, 0.5840584, 0.6933874, 0.7474869, 0.8722767, 0.9640924,
      1.0720271, 1.1024119, 1.1481851, 1.1753751, 1.2010460, 1.2520836, 1.3219817, 1.3748426, 1.6804600, 1.7541798,
      1.7597484, 1.8995670, 2.9941669, 3.1084058, 3.4349097, 3.6524039;

  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR(excitationEnergy_reference(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
