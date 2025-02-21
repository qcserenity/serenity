/**
 * @file NTOCalculator_test.cpp
 *
 * @date 29 Jan. 2021
 * @author Johannes Toelle
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
#include "postHF/LRSCF/Analysis/NTOCalculator.h"
#include "analysis/populationAnalysis/LoewdinPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h"
#include "data/ElectronicStructure.h"
#include "integrals/OneElectronIntegralController.h"
#include "misc/SystemSplittingTools.h"
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
 * @class NTOCalculatorTest
 * @brief Sets everything up for the tests of NTOCalculator.h/.cpp .
 */
class NTOCalculatorTest : public ::testing::Test {
 protected:
  NTOCalculatorTest() {
  }

  virtual ~NTOCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(NTOCalculatorTest, Super_NTOR) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(a, lrscf.settings);
  lrscf.settings.nEigen = 3;
  lrscf.run();
  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  // Values from Serenity 02.08.2019 // Cross checker with Turbomole 7.2.1 with RI and default Grid settings
  eigenvalues << 0.480578, 0.773088, 1.177130;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }
  // NTO calculator
  NTOCalculator<Options::SCF_MODES::RESTRICTED> ntos({a}, {}, 0.1);
  auto occEig = ntos.getOccEigenvalues(0);
  auto virtEig = ntos.getVirtEigenvalues(0);
  Eigen::VectorXd occEigRef(1);
  occEigRef << 0.909304006;
  Eigen::VectorXd virtEigRef(Eigen::VectorXd::Zero(9));
  virtEigRef[8] = 0.909304006;
  for_spin(occEig, virtEig) {
    for (unsigned int i = 0; i < virtEigRef.size(); ++i) {
      EXPECT_NEAR(virtEigRef(i), virtEig_spin(i), 1.0e-6);
    }
    for (unsigned int i = 0; i < occEigRef.size(); ++i) {
      EXPECT_NEAR(occEigRef(i), occEig_spin(i), 1.0e-6);
    }
  };
  auto occEig2 = ntos.getOccEigenvalues(1);
  auto virtEig2 = ntos.getVirtEigenvalues(1);
  Eigen::VectorXd occEigRef2(1);
  occEigRef2 << 8.8755256549e-01;
  Eigen::VectorXd virtEigRef2(Eigen::VectorXd::Zero(9));
  virtEigRef2[8] = 8.8755256549e-01;
  for_spin(occEig2, virtEig2) {
    for (unsigned int i = 0; i < virtEigRef2.size(); ++i) {
      EXPECT_NEAR(virtEigRef2(i), virtEig2_spin(i), 1.0e-6);
    }
    for (unsigned int i = 0; i < occEigRef2.size(); ++i) {
      EXPECT_NEAR(occEigRef2(i), occEig2_spin(i), 1.0e-6);
    }
  };

  std::remove("TestSystem_H2_def2_SVP_ACTIVE_FDE/NTOS/NTO1/README");
  std::remove("TestSystem_H2_def2_SVP_ACTIVE_FDE/NTOS/NTO1");
  std::remove("TestSystem_H2_def2_SVP_ACTIVE_FDE/NTOS/NTO2/README");
  std::remove("TestSystem_H2_def2_SVP_ACTIVE_FDE/NTOS/NTO2");
  std::remove("TestSystem_H2_def2_SVP_ACTIVE_FDE/NTOS");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(NTOCalculatorTest, SubCoupled_NTOR) {
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
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfA.settings.grid.accuracy = 7;
  lrscfA.settings.grid.smallGridAccuracy = 7;
  lrscfA.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.method = Options::LR_METHOD::TDA;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfB.settings.grid.accuracy = 7;
  lrscfB.settings.grid.smallGridAccuracy = 7;
  lrscfB.run();
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(act, lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(env, lrscfAB.settings);
  lrscfAB.settings.method = Options::LR_METHOD::TDA;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  lrscfAB.settings.grid.accuracy = 7;
  lrscfAB.settings.grid.smallGridAccuracy = 7;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.run();
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
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
  // NTO calculator
  NTOCalculator<Options::SCF_MODES::RESTRICTED> ntos({act, env}, {}, 0.1);
  auto occEig = ntos.getOccEigenvalues(0);
  auto virtEig = ntos.getVirtEigenvalues(0);
  Eigen::VectorXd occEigRef(2);
  occEigRef << 8.6676e-02, 0.91332431799220271;
  Eigen::VectorXd virtEigRef(12);
  virtEigRef << -8.0510e-19, -1.2906e-19, -6.2459e-21, 2.3264e-23, 3.7653e-21, 5.3121e-20, 9.6814e-20, 2.6685e-19,
      2.3357e-18, 2.4128e-17, 8.6676e-02, 0.91332430697171263;
  for_spin(occEig, virtEig) {
    for (unsigned int i = 0; i < virtEigRef.size(); ++i) {
      EXPECT_NEAR(virtEigRef(i), virtEig_spin(i), 1.0e-6);
    }
    for (unsigned int i = 0; i < occEigRef.size(); ++i) {
      EXPECT_NEAR(occEigRef(i), occEig_spin(i), 1.0e-6);
    }
  };
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/NTOS/NTO1/README");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/NTOS/NTO1");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/NTOS");
  std::remove("TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/NTOS/NTO1/README");
  std::remove("TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/NTOS/NTO1");
  std::remove("TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/NTOS");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(NTOCalculatorTest, Super_NTOU) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.grid.accuracy = 5;
  settings.grid.smallGridAccuracy = 3;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<Options::SCF_MODES::UNRESTRICTED>(system);
  scf.run();
  // TDDFT
  LRSCFTask<UNRESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<UNRESTRICTED>>(system, lrscf.settings);
  lrscf.settings.nEigen = 3;
  lrscf.settings.grid.accuracy = 5;
  lrscf.settings.grid.smallGridAccuracy = 5;
  lrscf.run();
  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  // NTO calculator
  NTOCalculator<Options::SCF_MODES::UNRESTRICTED> ntos({system}, {}, 0.1);
  auto occEig = ntos.getOccEigenvalues(0);
  auto virtEig = ntos.getVirtEigenvalues(0);
  Eigen::VectorXd occEigRef(2);
  occEigRef << .98967189000694378, 2.1462e-02;
  Eigen::VectorXd virtEigRef(2);
  virtEigRef << 0.98967189000694333, 2.1462e-02;
  unsigned counter = 0;
  for_spin(occEig, virtEig) {
    EXPECT_NEAR(virtEigRef(counter), virtEig_spin(virtEig_spin.size() - 1), 1.0e-6);
    EXPECT_NEAR(occEigRef(counter), occEig_spin(occEig_spin.size() - 1), 1.0e-6);
    counter++;
  };
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  std::remove("TestSystem_MethylRad_def2_SVP_PBE/NTOS/NTO1/README");
  std::remove("TestSystem_MethylRad_def2_SVP_PBE/NTOS/NTO1");
  std::remove("TestSystem_MethylRad_def2_SVP_PBE/NTOS");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(NTOCalculatorTest, Sub_Holedensity) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE_LRSCF);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT_LRSCF);
  LRSCFTaskSettings lrscfsettings;
  std::shared_ptr<LRSCFController<RESTRICTED>> lrscfContrAct =
      std::make_shared<LRSCFController<RESTRICTED>>(act, lrscfsettings);
  std::shared_ptr<LRSCFController<RESTRICTED>> lrscfContrEnv =
      std::make_shared<LRSCFController<RESTRICTED>>(env, lrscfsettings);
  auto excitationsActCoupled = lrscfContrAct->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv->getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  for (unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsEnvCoupled)(iExc), 1e-6);
  }
  // NTO calculator
  NTOCalculator<Options::SCF_MODES::RESTRICTED> ntoCalculator({act, env}, {}, 0.1);

  Settings settings = act->getSettings();
  settings.name = "TestSystem_H2Dimer_Def2_TZVP";
  auto supersys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP, settings);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({supersys});
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.nEigen = 20;
  lrscf.settings.conv = 1e-7;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.run();
  auto lrscfSuper = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(supersys, lrscf.settings);
  auto excitationsSuper = lrscfSuper->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  for (unsigned int iExc = 0; iExc < (*excitationsSuper).size(); iExc++) {
    EXPECT_NEAR((*excitationsActCoupled)(iExc), (*excitationsSuper)(iExc), 1e-6);
  }
  NTOCalculator<Options::SCF_MODES::RESTRICTED> superntos({supersys}, {}, 0.1);
  const MatrixInBasis<RESTRICTED>& superTrans = superntos.getTransitionDensity(0, 0);
  const MatrixInBasis<RESTRICTED>& superHole = superntos.getHoleDensity(0, 0);
  const MatrixInBasis<RESTRICTED>& superPart = superntos.getParticleDensity(0, 0);

  auto superBasisController = superTrans.getBasisController();

  MatrixInBasis<RESTRICTED> totalTrans =
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getTransitionDensity(0, 0), superBasisController) +
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getTransitionDensity(0, 1), superBasisController);
  MatrixInBasis<RESTRICTED> totalHole =
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getHoleDensity(0, 0), superBasisController) +
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getHoleDensity(0, 1), superBasisController);
  totalHole += SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getHoleDensityCorrection(0, 0),
                                                                           superBasisController) +
               SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getHoleDensityCorrection(0, 1),
                                                                           superBasisController);
  MatrixInBasis<RESTRICTED> totalPart =
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getParticleDensity(0, 0), superBasisController) +
      SystemSplittingTools<RESTRICTED>::projectMatrixIntoNewBasis(ntoCalculator.getParticleDensity(0, 1), superBasisController);

  const Eigen::MatrixXd transDiff = -1.0 * superTrans - totalTrans;
  const Eigen::MatrixXd holeDiff = superHole - totalHole;
  const Eigen::MatrixXd partDiff = superPart - totalPart;

  EXPECT_NEAR(transDiff.cwiseAbs().maxCoeff(), 0.0, 1e-6);
  EXPECT_NEAR(holeDiff.cwiseAbs().maxCoeff(), 0.0, 1e-6);
  EXPECT_NEAR(partDiff.cwiseAbs().maxCoeff(), 0.0, 1e-6);
  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/NTOS/NTO1/README");
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/NTOS/NTO1");
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/NTOS");
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/TestSystem_H2Dimer_Def2_TZVP/NTOS/NTO1/README");
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/TestSystem_H2Dimer_Def2_TZVP/NTOS/NTO1");
  std::remove("TestSystem_H2_MINBAS_ACTIVE_LRSCF/TestSystem_H2Dimer_Def2_TZVP/NTOS");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersys);
  std::remove("TestSystem_H2_MINBAS_ENVIRONMENT_LRSCF/NTOS/NTO1/README");
  std::remove("TestSystem_H2_MINBAS_ENVIRONMENT_LRSCF/NTOS/NTO1");
  std::remove("TestSystem_H2_MINBAS_ENVIRONMENT_LRSCF/NTOS");
}

} /* namespace Serenity */
