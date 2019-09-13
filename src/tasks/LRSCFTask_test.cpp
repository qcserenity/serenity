/**
 * @file LRSCFTask_test.cpp
 *
 * @date May 31, 2018
 * @author J. Toelle
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

/* Include Class Header*/
#include "tasks/LRSCFTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "tasks/FDETask.h"
#include "data/ElectronicStructure.h"
#include "tasks/LocalizationTask.h"
/* Include Serenity Internal Headers */
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "io/HDF5.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFTask.h/.cpp .
 */
class LRSCFTaskTest : public ::testing::Test {
 protected:
  LRSCFTaskTest() {
  }

  virtual ~LRSCFTaskTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(LRSCFTaskTest,TDDFT_LDA) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
 LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //Values from Serenity 02.08.2019 // Cross checker with Turbomole 7.2.1 with RI and default Grid settings
  eigenvalues <<   0.480578, 0.773088,1.177130;


  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  std::remove((a->getSettings().path+"H2_def2_SVP_ACTIVE_FDE.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();

}

TEST_F(LRSCFTaskTest,TDA_DFT) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.tda = true;
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //eigenvalues << 0.47220539945556966,
  //0.63852395474641876,
  //0.92676394912229387;

  eigenvalues<<0.47373563,0.63861229,0.92667761;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_PBE_NORI.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();

}
TEST_F(LRSCFTaskTest,CIS) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.tda = true;
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //These are the same values in ORCA 4_1_0
  eigenvalues<<0.50213139,0.64931589,0.94004596;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_HF.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(LRSCFTaskTest,CIS_UNRESTRICTED) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.tda = true;
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //These are the same values in ORCA 4_1_0
  eigenvalues<<0.36881612,0.50213139,0.52525269;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_HF_UNRESTRICTED.iso_lrscf.unres.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

//TDA Coupled + excldue Projection + BP86 + levelshift
TEST_F(LRSCFTaskTest,sTDA_Level) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
	
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
	task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task.run();
  //From Serenity 02.08.2019
	EXPECT_NEAR(-1.8155753808,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  auto task2 = TDEmbeddingTask<RESTRICTED>(env,act);
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.useEnvSys = true;
  task2.settings.enforceCharges = true;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task2.run();
  //From Serenity 02.08.2019
	EXPECT_NEAR(-1.8155753808,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);



  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act},{env});
  lrscfA.settings.tda = true;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env},{act});
  lrscfB.settings.tda = true;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act,env},{});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(env,lrscfAB.settings);
  
  lrscfAB.settings.tda = true;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.run();
  
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<<0.408893,
                              0.574051,
                              0.758401,
                              0.944049,
                              1.102164,
                              1.179569,
                              1.318655,
                              1.372784,
                              1.722485,
                              1.875651,
                              3.086483,
                              3.627045;

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.energies.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.orbs.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.dmat.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.settings").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.xyz").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.energies.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.orbs.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.dmat.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/out").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdec_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdec_lrscf.res.h5").c_str());

	SystemController__TEST_SUPPLY::cleanUp();
}

//TDA Coupled + excldue Projection + BP86 + Huzinaga
TEST_F(LRSCFTaskTest,sTDA_Huz) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
	
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
	task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task.run();

	EXPECT_NEAR(-1.8155753808,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  auto task2 = TDEmbeddingTask<RESTRICTED>(env,act);
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  task2.settings.useEnvSys = true;
  task2.settings.enforceCharges = true;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task2.run();

	EXPECT_NEAR(-1.8155753808,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act},{env});
  lrscfA.settings.tda = true;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env},{act});
  lrscfB.settings.tda = true;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act,env},{});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(env,lrscfAB.settings);
  
  lrscfAB.settings.tda = true;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  lrscfAB.run();
  
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<<0.408893,
                              0.574051,
                              0.758401,
                              0.944049,
                              1.102164,
                              1.179569,
                              1.318655,
                              1.372784,
                              1.722485,
                              1.875651,
                              3.086483,
                              3.627045;

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.energies.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.orbs.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.dmat.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.settings").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.xyz").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.energies.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.orbs.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.dmat.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/out").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdec_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdec_lrscf.res.h5").c_str());

	SystemController__TEST_SUPPLY::cleanUp();
}

//TDA Coupled + excldue Projection + BP86 + Hoffmann
TEST_F(LRSCFTaskTest,sTDA_Hof) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_BP86, true);
	
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
	task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task.run();

	EXPECT_NEAR(-1.8155753808,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  auto task2 = TDEmbeddingTask<RESTRICTED>(env,act);
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  task2.settings.useEnvSys = true;
  task2.settings.enforceCharges = true;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task2.run();

	EXPECT_NEAR(-1.8155753808,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act},{env});
  lrscfA.settings.tda = true;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env},{act});
  lrscfB.settings.tda = true;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act,env},{});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(env,lrscfAB.settings);
  
  lrscfAB.settings.tda = true;
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::HOFFMANN;
  lrscfAB.run();
  
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<<0.408893,
                              0.574051,
                              0.758401,
                              0.944049,
                              1.102164,
                              1.179569,
                              1.318655,
                              1.372784,
                              1.722485,
                              1.875651,
                              3.086483,
                              3.627045;

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.energies.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.orbs.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.dmat.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.settings").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.xyz").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.energies.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.orbs.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.dmat.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86/out").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86+TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_BP86.fdec_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdeu_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_BP86.fdec_lrscf.res.h5").c_str());

	SystemController__TEST_SUPPLY::cleanUp();
}


//TDDFT + Hybrid + NORI
TEST_F(LRSCFTaskTest,TDDFT_CAMB3LYP_NORI) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  eigenvalues<<0.47523804,0.63733218,0.92198322;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_CAMB3LYP.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}


//Everything with RI!

//TDHF
TEST_F(LRSCFTaskTest, TDHF){
    auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //These are values from ORCA 4_1_0
  eigenvalues<<0.497175,0.645541,0.931557;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_HF.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

//TDHF unrestricted
TEST_F(LRSCFTaskTest,TDHF_UNRESTRICTED) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.nEigen = 3;
  lrscf.run();
  
  auto excitations = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);
  
  //These are values from ORCA 4_1_0
  eigenvalues<<0.352212,0.497175,0.519376;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  std::remove((a->getSettings().path+"TestSystem_H2_DEF2_TZVP_HF_UNRESTRICTED.iso_lrscf.unres.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

//TDDFT unrestricted
TEST_F(LRSCFTaskTest,TDDFT_UNRESTRICTED) {
  auto a =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);
  
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto settings  = a->getSettings();
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED> >(a,lrscf.settings);
  
  lrscf.settings.nEigen = 5;
  lrscf.run();
  
  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(5);
  
  //These are values from Serenity 02.08.2019
  eigenvalues << 0.251192,
                0.254695,
                0.254849, 
                0.311227, 
                0.311302; 

  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-5);
  }

  std::remove((a->getSettings().path+"MethylRad_def2-SVP_PBE.iso_lrscf.unres.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

//Embedding mixed unrestricted restricted

//FDEc TDDFT + excldue Projection + B3LYP + levelshift
TEST_F(LRSCFTaskTest,sTDDFT_Level) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_B3LYP, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE_B3LYP, true);
	
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
	task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::B3LYP;
	task.run();

	EXPECT_NEAR(-1.8049132138,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  auto task2 = FDETask<RESTRICTED>(env,{act});
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::B3LYP;
	task2.run();

	EXPECT_NEAR(-1.8049132138,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act},{env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::B3LYP;
  lrscfA.settings.superSystemGrid = true;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env},{act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::B3LYP;
  lrscfA.settings.superSystemGrid = true;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act,env},{});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(env,lrscfAB.settings);
  
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::B3LYP;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfA.settings.superSystemGrid = true;
  lrscfAB.run();
  
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<<0.41363795, 
                              0.58140045,
                              0.75028406, 
                              0.95883377,
                              1.12931793, 
                              1.19638106, 
                              1.32308888, 
                              1.36862945, 
                              1.74765800, 
                              1.88834378, 
                              3.11266512, 
                              3.63866384;

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.energies.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.orbs.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.dmat.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.settings").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.xyz").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.energies.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.orbs.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.dmat.res.h5").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP/out").c_str());
	std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP+TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP").c_str());

  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.fdeu_lrscf.res.h5").c_str());
  std::remove((act->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE_B3LYP.fdec_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.fdeu_lrscf.res.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ENVIRONMENT_FDE_B3LYP.fdec_lrscf.res.h5").c_str());

	SystemController__TEST_SUPPLY::cleanUp();
}

//FDEc TDDFT + excldue Projection + PB86 + RI + levelshift + FUllFDEc
TEST_F(LRSCFTaskTest,sTDDFT_Level_RI) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_BP86, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_BP86, true);
  auto task = TDEmbeddingTask<RESTRICTED>(act,env);
	task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task.run();

	EXPECT_NEAR(-1.7736691992,act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  auto task2 = FDETask<RESTRICTED>(env,{act});
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
	task2.run();

	EXPECT_NEAR(-1.7736691992,env->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(),1e-6);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act},{env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfA.settings.superSystemGrid = true;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env},{act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfB.settings.superSystemGrid = true;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act,env},{});

  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(env,lrscfAB.settings);
  
  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.superSystemGrid = true;
  lrscfAB.settings.fullFDEc = true; 
  lrscfAB.run();
  
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation Serenity 05.08.2019
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<<0.338643,
                              0.474105,
                              0.590334,
                              0.855151,
                              0.897988,
                              0.968223,
                              1.054612,
                              1.113219,
                              1.116225,
                              1.405253,
                              1.509148,
                              1.767480;


  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }
	SystemController__TEST_SUPPLY::cleanUp();
}

//FDEc TDDFT + excldue Projection + PBE0 + RI + levelshift + FUllFDEc
TEST_F(LRSCFTaskTest,sTDDFT_Level_RI_Unrestricted) {
  auto act =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE_PBE0_UNRESTRICTED, true);
	auto env =
			SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ENVIRONMENT_FDE_PBE0_UNRESTRICTED, true);
  
  auto task = TDEmbeddingTask<UNRESTRICTED>(act,env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
	task.run();
  // This threshold is set to 1e-5 because the comparison with the supersystem energy is difficult
  // In the supersystem GS calculation no RI is used, in embedding case the coulomb term is however evaluated with RI
  // This also leads to larger deviation between supersystem and embedding TDDFT results later!!!
  // It was however checked that for NoRI both resutls are identical
	EXPECT_NEAR(-1.7561624580,act->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-5);

  auto task2 = FDETask<UNRESTRICTED>(env,{act});
	task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
	task2.run();

	EXPECT_NEAR(-1.7561624580,env->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(),1e-5);

  LRSCFTask<UNRESTRICTED> lrscfA({act},{env});
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.nEigen = 6;
  lrscfA.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
  lrscfA.settings.superSystemGrid = true;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfB({env},{act});
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.nEigen = 6;
  lrscfB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
  lrscfB.settings.superSystemGrid = true;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfAB({act,env},{});

  lrscfAB.settings.excludeProjection = true;
  lrscfAB.settings.nEigen = 12;
  lrscfAB.settings.embedding.naddXCFunc = Options::XCFUNCTIONALS::PBE0;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.superSystemGrid = true;
  lrscfAB.settings.fullFDEc = true; 
  lrscfAB.run();
  
  auto lrscfContrAct = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED> >(act,lrscfAB.settings);
  auto lrscfContrEnv = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED> >(env,lrscfAB.settings);
  auto excitationsActUncoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsActCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);
  auto excitationsEnvUnoupled = lrscfContrAct -> getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  auto excitationsEnvCoupled = lrscfContrEnv -> getExcitationEnergies(Options::LRSCF_TYPE::COUPLED);

  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR((*excitationsActCoupled)(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

  //From Supersystem calculation Serenity 06.08.2019
  Eigen::VectorXd excitationEnergy_reference(12);
  excitationEnergy_reference<< 0.282916,
                                0.352153,
                                0.420314,
                                0.489925,
                                0.549325,
                                0.592604,
                                0.826059,
                                0.882117,
                                0.985307,
                                1.007757,
                                1.101853,
                                1.132371;


  for(unsigned int iExc = 0; iExc < (*excitationsActCoupled).size(); iExc ++){
    EXPECT_NEAR(excitationEnergy_reference(iExc),(*excitationsEnvCoupled)(iExc),1e-6);
  }

	SystemController__TEST_SUPPLY::cleanUp();
}

// ToDo FDEc non-additive kinetic

//local Mos
TEST_F(LRSCFTaskTest,TDDFT_LOCAL_MOs) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE_BP86, true);

  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(a);

  //First in canonical MOs and then in local MOs  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscf.settings);
  lrscf.settings.nEigen = 3;
  lrscf.run();
  auto excitationsCanon = lrscfContr -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);

  //Localization
  auto task = LocalizationTask(a);
	task.run();

  //First in canonical MOs and then in local MOs  
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfNonCanon({a});
  auto lrscfContrNonCanon = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED> >(a,lrscfNonCanon.settings);
  lrscfNonCanon.settings.nEigen = 3;
  lrscfNonCanon.settings.localMO = true;
  lrscfNonCanon.run();
  auto excitationsNonCanon = lrscfContrNonCanon -> getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);

 for(unsigned int iExc = 0; iExc < (*excitationsNonCanon).size(); iExc ++){
    EXPECT_NEAR((*excitationsNonCanon)(iExc),(*excitationsCanon)(iExc),1e-7);
  }
  std::remove((a->getSettings().path+"TestSystem_He2_6_31Gs_BP86.iso_lrscf.res.h5").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}


}
