/**
 * @file LRSCFTaskBSE_test.cpp
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
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GWTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFRITask.h/.cpp .
 */
class LRSCFTaskBSETest : public ::testing::Test {
 protected:
  LRSCFTaskBSETest() {
  }
  virtual ~LRSCFTaskBSETest() = default;
  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(LRSCFTaskBSETest, RPAScreeningTDA_R) {
  auto a_noRI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  Settings settings = a_noRI->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  auto a_RI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, settings);
  // SCF
  ScfTask<RESTRICTED> scf_RI(a_RI);
  scf_RI.run();
  // LRSCF
  std::vector<std::shared_ptr<SystemController>> active_RI;
  active_RI.push_back(a_RI);
  LRSCFTask<RESTRICTED> lrscf_RI(active_RI);
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(a_RI, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 3;
  lrscf_RI.settings.method = Options::LR_METHOD::TDA;
  lrscf_RI.run();
  auto excitations_RI = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);

  // Values from old BSE code
  Eigen::VectorXd ref(3);
  ref << 0.185510, 0.344501, 0.583874;

  for (unsigned int i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(ref(i), (*excitations_RI)(i), 1.0e-6);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_noRI);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_RI);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskBSETest, RPAScreeningTDDFT_R) {
  auto a_noRI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  Settings settings = a_noRI->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  auto a_RI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, settings);
  // SCF
  ScfTask<RESTRICTED> scf_RI(a_RI);
  scf_RI.run();
  // LRSCF
  std::vector<std::shared_ptr<SystemController>> active_RI;
  active_RI.push_back(a_RI);
  LRSCFTask<RESTRICTED> lrscf_RI(active_RI);
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(a_RI, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 3;
  lrscf_RI.run();
  auto excitations_RI = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);

  // Values from old BSE code
  Eigen::VectorXd ref(3);
  ref << 0.145630, 0.331067, 0.552855;
  for (unsigned int i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(ref(i), (*excitations_RI)(i), 1.0e-6);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_noRI);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_RI);
  SystemController__TEST_SUPPLY::cleanUp();
}

/* TEST_F(LRSCFTaskBSETest, RPAScreeningTDDFT_R_NAF) {
  auto a_noRI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  Settings settings = a_noRI->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.densityFitting = Options::DENS_FITS::RI;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  auto a_RI = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, settings);
  // SCF
  ScfTask<RESTRICTED> scf_RI(a_RI);
  scf_RI.run();
  // LRSCF
  std::vector<std::shared_ptr<SystemController>> active_RI;
  active_RI.push_back(a_RI);
  LRSCFTask<RESTRICTED> lrscf_RI(active_RI);
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(a_RI, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = true;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 3;
  lrscf_RI.settings.nafThresh = 1e-2;
  lrscf_RI.run();
  auto excitations_RI = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);

  // Values from old BSE code
  Eigen::VectorXd ref(3);
  ref << 0.145630, 0.331067, 0.552855;
  for (unsigned int i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(ref(i), (*excitations_RI)(i), 1.0e-6);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_noRI);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a_RI);
  SystemController__TEST_SUPPLY::cleanUp();
} */

TEST_F(LRSCFTaskBSETest, BSE_TDA) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  std::vector<std::shared_ptr<SystemController>> active;

  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(a);
  scf.run();
  // Turbomole 7.4.1 ridft/s-vwn/m3
  EXPECT_NEAR(-1.13176185536, a->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);

  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({a});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 100;
  gw.settings.nOcc = 100;
  gw.settings.linearized = false;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.00001;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/CD/eta 0.001
  Eigen::VectorXd correlationEnergyTurbo_eV(10);
  correlationEnergyTurbo_eV << -0.188, -0.238, -0.561, -2.910, -1.599, -1.599, -2.164, -0.061, -0.061, -0.184;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 1.5e-3);
    }
  };
  LRSCFTask<RESTRICTED> lrscf_RI({a});
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(a, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDA;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 10;
  lrscf_RI.run();
  auto excitationEnergies = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  // Turbomole 7.4.1 BSE/ri approximation/all states
  Eigen::VectorXd bseExcitationEnergyTurbo_eV(9);
  bseExcitationEnergyTurbo_eV << 14.77514, 22.16088, 31.58836, 43.96886, 43.96886, 59.99772, 64.68772, 64.68772, 106.32841;

  for (unsigned int i = 0; i < bseExcitationEnergyTurbo_eV.size(); i++) {
    EXPECT_NEAR((*excitationEnergies)(i)*HARTREE_TO_EV, bseExcitationEnergyTurbo_eV(i), 1.7e-3);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskBSETest, BSE_TDDFT) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  std::vector<std::shared_ptr<SystemController>> active;

  // Perform SCF
  auto scf = ScfTask<RESTRICTED>(a);
  scf.run();
  // Turbomole 7.4.1 ridft/s-vwn/m3
  EXPECT_NEAR(-1.13176185536, a->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(), 1e-5);

  // Perform CD GW calculation
  auto gw = GWTask<RESTRICTED>({a});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 100;
  gw.settings.nOcc = 100;
  gw.settings.linearized = false;
  // set to zero because closer to the turbomole results
  gw.settings.eta = 0.00001;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  // Turbomole 7.4.1 gw/ri approximation/CD/eta 0.001
  Eigen::VectorXd correlationEnergyTurbo_eV(10);
  correlationEnergyTurbo_eV << -0.188, -0.238, -0.561, -2.910, -1.599, -1.599, -2.164, -0.061, -0.061, -0.184;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlationEnergyTurbo_eV.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlationEnergyTurbo_eV(i), 1.5e-3);
    }
  };
  LRSCFTask<RESTRICTED> lrscf_RI({a});
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(a, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 10;
  lrscf_RI.run();
  auto excitationEnergies = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  // Turbomole 7.4.1 BSE/ri approximation/all states
  Eigen::VectorXd bseExcitationEnergyTurbo_eV(9);
  bseExcitationEnergyTurbo_eV << 14.38243, 21.94162, 31.12418, 43.78195, 43.78195, 59.92243, 64.66504, 64.66504, 106.29828;

  for (unsigned int i = 0; i < bseExcitationEnergyTurbo_eV.size(); i++) {
    EXPECT_NEAR((*excitationEnergies)(i)*HARTREE_TO_EV, bseExcitationEnergyTurbo_eV(i), 1.7e-3);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskBSETest, BSE_TDDFT_U) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  // Perform SCF
  auto scf = ScfTask<UNRESTRICTED>(system);
  scf.run();
  // Perform CD GW calculation
  auto gw = GWTask<UNRESTRICTED>({system});
  gw.settings.gwtype = Options::GWALGORITHM::CD;
  gw.settings.nVirt = 100;
  gw.settings.nOcc = 100;
  gw.settings.linearized = false;
  gw.settings.eta = 0.0;
  gw.settings.integrationPoints = 1000;
  gw.run();
  LRSCFTask<UNRESTRICTED> lrscf_RI({system});
  auto lrscfContr_RI = std::make_shared<LRSCFController<UNRESTRICTED>>(system, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 10;
  lrscf_RI.run();
  auto excitationEnergies = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  // Turbomole 7.4.1 BSE/ri approximation/all states
  Eigen::VectorXd bseExcitationEnergyTurbo_eV(10);
  bseExcitationEnergyTurbo_eV << 7.64654, 8.04465, 8.04959, 9.41092, 9.41305, 12.06516, 12.07026, 13.02549, 13.59712, 13.60258;

  for (unsigned int i = 0; i < bseExcitationEnergyTurbo_eV.size(); i++) {
    EXPECT_NEAR((*excitationEnergies)(i)*HARTREE_TO_EV, bseExcitationEnergyTurbo_eV(i), 2.9e-3);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskBSETest, embBSE) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);

  Settings settings1 = act->getSettings();
  settings1.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings1.basis.label = "DEF2-TZVP";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP_B2PLYP, settings1);

  Settings settings2 = env->getSettings();
  settings2.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings2.basis.label = "DEF2-TZVP";
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, settings2);

  auto task = FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>({act, env}, {});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  task.run();

  auto gw = GWTask<Options::SCF_MODES::RESTRICTED>({act}, {env});
  gw.settings.eta = 0.001;
  gw.settings.diis = true;
  gw.settings.nVirt = 1;
  gw.settings.nOcc = 1;
  gw.settings.qpiterations = 10;
  gw.settings.ConvergenceThreshold = 1e-6;
  gw.settings.integrationPoints = 128;
  gw.settings.subsystemAuxillaryBasisOnly = true;
  gw.run();
  auto correlationEnergy = gw._correlationEnergies;
  Eigen::VectorXd correlation2(6);
  // Serenity 19.05.2022 (from CD)
  correlation2 << 0.0, 0.0, 0.0, 0.0, 1.819371555985553, -0.69471338695164397;
  for_spin(correlationEnergy) {
    for (unsigned int i = 0; i < correlation2.size(); i++) {
      EXPECT_NEAR(correlationEnergy_spin(i) * HARTREE_TO_EV, correlation2(i), 2e-6);
    }
  };
  LRSCFTask<RESTRICTED> lrscf_RI({act}, {env});
  auto lrscfContr_RI = std::make_shared<LRSCFController<RESTRICTED>>(act, lrscf_RI.settings);
  lrscf_RI.settings.method = Options::LR_METHOD::TDDFT;
  lrscf_RI.settings.densFitK = Options::DENS_FITS::RI;
  lrscf_RI.settings.rpaScreening = true;
  lrscf_RI.settings.nEigen = 10;
  lrscf_RI.run();
  auto excitationEnergies = lrscfContr_RI->getExcitationEnergies(Options::LRSCF_TYPE::UNCOUPLED);
  // Turbomole 7.4.1 BSE/ri approximation/all states
  Eigen::VectorXd bseExcitationEnergy_eV(10);
  bseExcitationEnergy_eV << 6.85445, 6.86877, 7.39249, 7.57497, 10.27804, 10.89427, 13.29167, 13.55584, 14.18332, 14.23498;
  for (unsigned int i = 0; i < bseExcitationEnergy_eV.size(); i++) {
    EXPECT_NEAR((*excitationEnergies)(i)*HARTREE_TO_EV, bseExcitationEnergy_eV(i), 1.0e-5);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
