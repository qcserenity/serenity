/**
 * @file LRSCFTaskTDDFT_test.cpp
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
#include "math/FloatMaths.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/LocalizationTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <tuple>

namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFTask.h/.cpp .
 */
class LRSCFTaskTDDFTTest : public ::testing::Test {
 protected:
  LRSCFTaskTDDFTTest() {
  }

  virtual ~LRSCFTaskTDDFTTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(LRSCFTaskTDDFTTest, TDDFT) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.nEigen = 3;
  lrscf.settings.frequencies = {1, 2, 3, 4, 5};
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // Values from Serenity 02.08.2019 // Cross checker with Turbomole 7.2.1 with RI and default Grid settings
  eigenvalues << 0.480578, 0.773088, 1.177130;

  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDA) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.nEigen = 3;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // eigenvalues << 0.47220539945556966,
  // 0.63852395474641876,
  // 0.92676394912229387;

  eigenvalues << 0.47373563, 0.63861229, 0.92667761;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, CIS) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 3;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // These are the same values in ORCA 4_1_0
  eigenvalues << 0.50213139, 0.64931589, 0.94004596;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, CIS_UNRESTRICTED) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 3;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // These are the same values in ORCA 4_1_0
  eigenvalues << 0.36881612, 0.50213139, 0.52525269;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-7);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_CAMB3LYP_NORI) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.nEigen = 3;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  eigenvalues << 0.47523804, 0.63733218, 0.92198322;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDHF) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.nEigen = 3;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // These are values from ORCA 4_1_0
  eigenvalues << 0.497175, 0.645541, 0.931557;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDHF_UNRESTRICTED) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF_UNRESTRICTED, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.nEigen = 3;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(3);

  // These are values from ORCA 4_1_0
  eigenvalues << 0.352212, 0.497175, 0.519376;
  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_UNRESTRICTED) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(sys, lrscf.settings);

  lrscf.settings.nEigen = 5;
  lrscf.run();

  auto excitations = lrscfContr->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED);
  Eigen::VectorXd eigenvalues(5);

  // These are values from Serenity 02.08.2019
  eigenvalues << 0.251192, 0.254695, 0.254849, 0.311227, 0.311302;

  for (unsigned int i = 0; i < (unsigned)lrscf.settings.nEigen; ++i) {
    EXPECT_NEAR(eigenvalues(i), (*excitations)(i), 1.0e-5);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, LMO_TDDFT) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);

  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  lrscf.settings.conv = 1E-6;
  lrscf.run();

  auto task = LocalizationTask(sys);
  task.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfNonCanon(active);
  lrscfNonCanon.settings.conv = 1E-6;
  lrscfNonCanon.run();

  double maxDiff = (lrscf.getTransitions() - lrscfNonCanon.getTransitions()).cwiseAbs().maxCoeff();
  EXPECT_LE(maxDiff, 1E-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, Response_Properties) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ, settings);

  // Run in undamped mode.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf1({sys});
  lrscf1.settings.nEigen = 3;
  lrscf1.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf1.settings.grid.blockAveThreshold = 1e-12;
  lrscf1.run();
  // Also run in damped mode.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf2({sys});
  lrscf2.settings.nEigen = 3;
  lrscf2.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf2.settings.grid.blockAveThreshold = 1e-12;
  lrscf2.settings.damping = 0.5;
  lrscf2.run();
  // Run in undamped mode (velocity gauge).
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf3({sys});
  lrscf3.settings.nEigen = 3;
  lrscf3.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf3.settings.grid.blockAveThreshold = 1e-12;
  lrscf3.settings.gauge = Options::GAUGE::VELOCITY;
  lrscf3.run();

  Eigen::MatrixXd excitations(lrscf1.settings.nEigen, 6);
  excitations.row(0) << 0.202326844, 0.000036960, 0.000032958, -0.000034701, -0.000032769, -0.000034701;
  excitations.row(1) << 0.211120996, 0.017366910, 0.017608909, -0.000184005, -0.000184698, -0.000183425;
  excitations.row(2) << 0.224485241, 0.002012985, 0.001853583, -0.000147352, -0.000141397, -0.000147352;

  Eigen::MatrixXd undamped_properties(lrscf1.settings.frequencies.size(), 5);
  undamped_properties.row(0) << 0.110247967, 32.586205728, -0.281570513, 0.000000000, 0.000000000;
  undamped_properties.row(1) << 0.146997290, 34.528141086, -0.491967329, 0.000000000, 0.000000000;
  undamped_properties.row(2) << 0.183746612, 38.513272661, -1.529629721, 0.000000000, 0.000000000;

  Eigen::MatrixXd damped_properties(lrscf2.settings.frequencies.size(), 5);
  damped_properties.row(0) << 0.110247967, 32.489845186, -0.258837132, 0.007232187, -0.023923107;
  damped_properties.row(1) << 0.146997290, 34.324869878, -0.413423700, 0.016916436, -0.092620820;
  damped_properties.row(2) << 0.183746612, 37.567563674, -0.752677381, 0.047185222, -0.618886578;

  Eigen::MatrixXd undamped_properties_velocity(lrscf3.settings.frequencies.size(), 5);
  undamped_properties_velocity.row(0) << 0.110247967, 31.9955162, -0.2642771, 0.000000000, 0.000000000;
  undamped_properties_velocity.row(1) << 0.146997290, 33.9248452, -0.4669986, 0.000000000, 0.000000000;
  undamped_properties_velocity.row(2) << 0.183746612, 37.8981875, -1.4782757, 0.000000000, 0.000000000;

  // Compare excitation energies and oscillator/rotator strengths.
  EXPECT_LE((excitations - lrscf1.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);
  EXPECT_LE((excitations - lrscf2.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);
  EXPECT_LE((excitations - lrscf3.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);

  // Compare undamped properties.
  for (unsigned iFreq = 0; iFreq < 3; ++iFreq) {
    double frequency = lrscf1.settings.frequencies[iFreq];
    double isoPol = 1. / 3. * std::get<1>(lrscf1.getProperties()[iFreq]).trace();
    double isoOR = -1. / 3. / frequency * std::get<2>(lrscf1.getProperties()[iFreq]).trace();
    double isoAbs = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf1.getProperties()[iFreq]).trace();
    double isoCD = 6.533 * frequency * std::get<4>(lrscf1.getProperties()[iFreq]).trace();
    EXPECT_NEAR(isoPol, undamped_properties(iFreq, 1), LOOSE_D);
    EXPECT_NEAR(isoOR, undamped_properties(iFreq, 2), LOOSE_D);
    EXPECT_NEAR(isoAbs, undamped_properties(iFreq, 3), LOOSE_D);
    EXPECT_NEAR(isoCD, undamped_properties(iFreq, 4), LOOSE_D);
  }
  // Compare undamped properties (velocity).
  for (unsigned iFreq = 0; iFreq < 3; ++iFreq) {
    double frequency = lrscf3.settings.frequencies[iFreq];
    double isoPol = 1. / 3. * std::get<1>(lrscf3.getProperties()[iFreq]).trace();
    double isoOR = -1. / 3. / frequency * std::get<2>(lrscf3.getProperties()[iFreq]).trace();
    double isoAbs = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf3.getProperties()[iFreq]).trace();
    double isoCD = 6.533 * frequency * std::get<4>(lrscf3.getProperties()[iFreq]).trace();
    EXPECT_NEAR(isoPol, undamped_properties_velocity(iFreq, 1), LOOSE_D);
    EXPECT_NEAR(isoOR, undamped_properties_velocity(iFreq, 2), LOOSE_D);
    EXPECT_NEAR(isoAbs, undamped_properties_velocity(iFreq, 3), LOOSE_D);
    EXPECT_NEAR(isoCD, undamped_properties_velocity(iFreq, 4), LOOSE_D);
  }
  // Compare damped properties.
  for (unsigned iFreq = 0; iFreq < 3; ++iFreq) {
    double frequency = lrscf2.settings.frequencies[iFreq];
    double isoPol = 1. / 3. * std::get<1>(lrscf2.getProperties()[iFreq]).trace();
    double isoOR = -1. / 3. / frequency * std::get<2>(lrscf2.getProperties()[iFreq]).trace();
    double isoAbs = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf2.getProperties()[iFreq]).trace();
    double isoCD = 6.533 * frequency * std::get<4>(lrscf2.getProperties()[iFreq]).trace();
    EXPECT_NEAR(isoPol, damped_properties(iFreq, 1), 1.5e-6);
    EXPECT_NEAR(isoOR, damped_properties(iFreq, 2), LOOSE_D);
    EXPECT_NEAR(isoAbs, damped_properties(iFreq, 3), LOOSE_D);
    EXPECT_NEAR(isoCD, damped_properties(iFreq, 4), LOOSE_D);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, Restart) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.basis.label = "DEF2-SVP";

  // Have to set this to zero so the test isn't trying to look
  // there for old vectors but in the current working directory instead.
  settings.load = "";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  unsigned smallEigen = 6;
  unsigned largeEigen = smallEigen * 2;

  // Initial run: sloppy convergence.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf1({sys});
  lrscf1.settings.nEigen = smallEigen;
  lrscf1.settings.conv = 1e-3;
  lrscf1.run();

  // Load the same number of vectors from disk and converge them.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf2({sys});
  lrscf2.settings.nEigen = smallEigen;
  lrscf2.settings.conv = 1e-6;
  lrscf2.settings.restart = true;
  lrscf2.run();

  // Start from scratch another time and compare.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf3({sys});
  lrscf3.settings.nEigen = smallEigen;
  lrscf3.settings.conv = 1e-6;
  lrscf3.run();

  EXPECT_LE((lrscf2.getTransitions() - lrscf3.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);

  // Load vectors from disk, but this time, they are not enough.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf4({sys});
  lrscf4.settings.nEigen = largeEigen;
  lrscf4.settings.conv = 1e-6;
  lrscf4.settings.restart = true;
  lrscf4.run();

  // Run from scratch and compare.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf5({sys});
  lrscf5.settings.nEigen = largeEigen;
  lrscf5.settings.conv = 1e-6;
  lrscf5.run();

  EXPECT_LE((lrscf4.getTransitions() - lrscf5.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);

  // Test if fully converged roots are still fully converged in 1 cycle.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf6({sys});
  lrscf6.settings.nEigen = largeEigen;
  lrscf6.settings.conv = 1e-6;
  lrscf6.settings.maxCycles = 1;
  lrscf6.settings.restart = true;
  lrscf6.run();

  EXPECT_LE((lrscf5.getTransitions() - lrscf6.getTransitions()).cwiseAbs().maxCoeff(), LOOSE_D);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, simplifiedTDA_Restricted) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.grimme = true;
  lrscfTDDFT.run();

  auto cont = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscfTDDFT.settings);
  auto excitations = *(cont->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED));

  Eigen::VectorXd referenceTDDFT(4);
  referenceTDDFT << 0.1446055, 0.2204230, 0.2427796, 0.2516147;

  for (unsigned i = 0; i < referenceTDDFT.size(); ++i) {
    EXPECT_LE(std::abs(referenceTDDFT[i] - excitations[i]), LOOSE_D);
  }

  // TDA
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDA({sys});
  lrscfTDA.settings.grimme = true;
  lrscfTDA.settings.method = Options::LR_METHOD::TDA;
  lrscfTDA.run();

  cont = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(sys, lrscfTDA.settings);
  excitations = *(cont->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED));

  Eigen::VectorXd referenceTDA(4);
  referenceTDA << 0.1447616, 0.2204958, 0.2428295, 0.2516237;

  for (unsigned i = 0; i < referenceTDA.size(); ++i) {
    EXPECT_LE(std::abs(referenceTDA[i] - excitations[i]), LOOSE_D);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, simplifiedTDA_Unrestricted) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.grimme = true;
  lrscfTDDFT.settings.nEigen = 8;
  lrscfTDDFT.run();

  auto cont = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(sys, lrscfTDDFT.settings);
  auto excitations = *(cont->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED));

  Eigen::VectorXd referenceTDDFT(8);
  referenceTDDFT << 0.1446055, 0.1446055, 0.2195534, 0.2204230, 0.2422150, 0.2427796, 0.2515966, 0.2516147;

  for (unsigned i = 0; i < referenceTDDFT.size(); ++i) {
    EXPECT_LE(std::abs(referenceTDDFT[i] - excitations[i]), LOOSE_D);
  }

  // TDA
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscfTDA({sys});
  lrscfTDA.settings.grimme = true;
  lrscfTDA.settings.method = Options::LR_METHOD::TDA;
  lrscfTDA.settings.nEigen = 8;
  lrscfTDA.run();

  cont = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(sys, lrscfTDA.settings);
  excitations = *(cont->getExcitationEnergies(Options::LRSCF_TYPE::ISOLATED));

  Eigen::VectorXd referenceTDA(8);
  referenceTDA << 0.1447616, 0.1447616, 0.2195566, 0.2204958, 0.2422169, 0.2428295, 0.2516039, 0.2516238;

  for (unsigned i = 0; i < referenceTDA.size(); ++i) {
    EXPECT_LE(std::abs(referenceTDA[i] - excitations[i]), LOOSE_D);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, DoubleHybridTDDFT) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.run();

  auto tddft = lrscfTDDFT.getTransitions();

  Eigen::MatrixXd refTDDFT(4, 6);
  refTDDFT.row(0) << 0.147861075, 0.000000000, 0.000000002, 0.000000169, 0.000000259, 0.000000093;
  refTDDFT.row(1) << 0.338190095, 0.001289817, 0.004528878, 0.000001611, 0.000001609, 0.000000859;
  refTDDFT.row(2) << 0.303410187, 0.149607420, 0.176382875, -0.000001683, -0.000001913, -0.000001762;
  refTDDFT.row(3) << 0.361563767, 0.076219956, 0.053208231, -0.000000174, -0.000000146, -0.000000174;
  EXPECT_LE((refTDDFT - tddft).cwiseAbs().maxCoeff(), 1e-6);

  // TDA
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDA({sys});
  lrscfTDA.settings.method = Options::LR_METHOD::TDA;
  lrscfTDA.run();
  auto tda = lrscfTDA.getTransitions();

  Eigen::MatrixXd refTDA(4, 6);
  refTDA.row(0) << 0.149730729, 0.000000000, 0.000000002, 0.000000185, 0.000000049, 0.000000022;
  refTDA.row(1) << 0.303456079, 0.163986142, 0.133651532, 0.000000665, 0.000000537, 0.000000595;
  refTDA.row(2) << 0.342005361, 0.001887603, 0.000000380, -0.000000723, -0.000000556, -0.000039189;
  refTDA.row(3) << 0.369321362, 0.044015294, 0.008240868, -0.000000159, 0.000000069, 0.000000159;
  EXPECT_LE((refTDA - tda).cwiseAbs().maxCoeff(), 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, DoubleHybridTDDFT_CustomFunctional) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                          BasicFunctionals::BASIC_FUNCTIONALS::C_LYP};
  settings.customFunc.mixingFactors = {0.47, 0.73};
  settings.customFunc.hfExchangeRatio = 0.53;
  settings.customFunc.hfCorrelRatio = 0.27;
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                     BasicFunctionals::BASIC_FUNCTIONALS::C_LYP};
  lrscfTDDFT.settings.customFunc.mixingFactors = {0.47, 0.73};
  lrscfTDDFT.settings.customFunc.hfExchangeRatio = 0.53;
  lrscfTDDFT.settings.customFunc.hfCorrelRatio = 0.27;
  lrscfTDDFT.settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfTDDFT.run();

  auto tddft = lrscfTDDFT.getTransitions();

  Eigen::MatrixXd refTDDFT(4, 6);
  refTDDFT.row(0) << 0.147861075, 0.000000000, 0.000000002, 0.000000169, 0.000000259, 0.000000093;
  refTDDFT.row(1) << 0.338190095, 0.001289817, 0.004528878, 0.000001611, 0.000001609, 0.000000859;
  refTDDFT.row(2) << 0.303410187, 0.149607420, 0.176382875, -0.000001683, -0.000001913, -0.000001762;
  refTDDFT.row(3) << 0.361563767, 0.076219956, 0.053208231, -0.000000174, -0.000000146, -0.000000174;
  EXPECT_LE((refTDDFT - tddft).cwiseAbs().maxCoeff(), 1e-6);

  // TDA
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDA({sys});
  lrscfTDA.settings.method = Options::LR_METHOD::TDA;
  lrscfTDA.settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                                   BasicFunctionals::BASIC_FUNCTIONALS::C_LYP};
  lrscfTDA.settings.customFunc.mixingFactors = {0.47, 0.73};
  lrscfTDA.settings.customFunc.hfExchangeRatio = 0.53;
  lrscfTDA.settings.customFunc.hfCorrelRatio = 0.27;
  lrscfTDA.settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  lrscfTDA.run();
  auto tda = lrscfTDA.getTransitions();

  Eigen::MatrixXd refTDA(4, 6);
  refTDA.row(0) << 0.149730729, 0.000000000, 0.000000002, 0.000000185, 0.000000049, 0.000000022;
  refTDA.row(1) << 0.303456079, 0.163986142, 0.133651532, 0.000000665, 0.000000537, 0.000000595;
  refTDA.row(2) << 0.342005361, 0.001887603, 0.000000380, -0.000000723, -0.000000556, -0.000039189;
  refTDA.row(3) << 0.369321362, 0.044015294, 0.008240868, -0.000000159, 0.000000069, 0.000000159;
  EXPECT_LE((refTDA - tda).cwiseAbs().maxCoeff(), 1e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_SmallLargeGrid) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.grid.smallGridAccuracy = 2;
  lrscfTDDFT.run();

  auto excitations = lrscfTDDFT.getTransitions();

  Eigen::MatrixXd reference(4, 6);
  reference.row(0) << 0.147312398818, 0.000000000100, 0.000000001919, 0.000000102957, 0.000000279110, 0.000000063801;
  reference.row(1) << 0.315462818892, 0.143124340930, 0.126260383138, -0.000000030699, -0.000000041334, -0.000000044008;
  reference.row(2) << 0.337621554713, 0.001438729439, 0.004959190605, -0.000000133653, -0.000000208838, -0.000000112485;
  reference.row(3) << 0.358724408102, 0.005999079321, 0.001096777113, -0.000000173967, -0.000000074580, -0.000000174423;

  EXPECT_LE((reference - excitations).cwiseAbs().maxCoeff(), LOOSE_D);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, FrozenCoreAndFrozenVirtual) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.frozenCore = true;
  lrscfTDDFT.settings.frozenVirtual = 20;
  lrscfTDDFT.run();

  auto tddft = lrscfTDDFT.getTransitions();

  Eigen::MatrixXd ref(4, 6);
  ref.row(0) << 0.144514342492, 0.000000001789, 0.000000002126, 0.000000070920, 0.000000001388, 0.000000001273;
  ref.row(1) << 0.281169485317, 0.121419972289, 0.100597633838, 0.000000005516, -0.000000002240, -0.000000002460;
  ref.row(2) << 0.334127018547, 0.002541983319, 0.006788195441, -0.000000030490, -0.000000023400, -0.000000014319;
  ref.row(3) << 0.342987679763, 0.028816932705, 0.044801409237, -0.000000027119, -0.000000033829, -0.000000027131;
  EXPECT_LE((ref - tddft).cwiseAbs().maxCoeff(), LOOSE_D);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TripletExcitations) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.label = "DEF2-SVP";

  // Test 1: LDA/TDA
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
  auto sys1 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  LRSCFTask<RESTRICTED> tda({sys1});
  tda.settings.triplet = true;
  tda.settings.method = Options::LR_METHOD::TDA;
  tda.run();

  Eigen::VectorXd tda_ref(4);
  tda_ref << 0.1150766, 0.2377234, 0.2624079, 0.2834016;
  EXPECT_LE((tda.getTransitions().col(0) - tda_ref).cwiseAbs().maxCoeff(), 5e-6);

  // Test 2: PBE/TDDFT
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  auto sys2 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  LRSCFTask<RESTRICTED> tddft({sys2});
  tddft.settings.triplet = true;
  tddft.run();

  Eigen::VectorXd tddft_ref(4);
  tddft_ref << 0.1143974, 0.2155707, 0.2570795, 0.2798955;
  EXPECT_LE((tddft.getTransitions().col(0) - tddft_ref).cwiseAbs().maxCoeff(), 5e-6);

  // Test 3: PBE0/RPA
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  auto sys3 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  LRSCFTask<RESTRICTED> rpa({sys3});
  rpa.settings.triplet = true;
  rpa.run();

  Eigen::VectorXd rpa_ref(4);
  rpa_ref << 0.1164686, 0.1945099, 0.2865134, 0.2883944;
  EXPECT_LE((rpa.getTransitions().col(0) - rpa_ref).cwiseAbs().maxCoeff(), 5e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_CoreOnly) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfTDDFT({sys});
  lrscfTDDFT.settings.coreOnly = true;
  lrscfTDDFT.run();

  auto excitations = lrscfTDDFT.getTransitions().col(0);

  Eigen::VectorXd referenceTDDFT(4);
  referenceTDDFT << 9.9070248, 9.9841180, 10.0111261, 10.0210513;

  EXPECT_LE((referenceTDDFT - excitations).cwiseAbs().maxCoeff(), LOOSE_D);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_SCFStab) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.basis.label = "DEF2-SVP";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> real({sys});
  real.settings.scfstab = Options::STABILITY_ANALYSIS::REAL;
  real.settings.densFitJ = Options::DENS_FITS::NONE;
  real.run();
  double result = real.getTransitions()(0, 0);
  // Turbomole 7.5.1 Feb 2023
  double reference = 0.1603296108583212;
  EXPECT_LE(std::abs(result - reference), 1e-4);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> nonreal({sys});
  nonreal.settings.scfstab = Options::STABILITY_ANALYSIS::NONREAL;
  nonreal.settings.densFitJ = Options::DENS_FITS::NONE;
  nonreal.run();
  result = nonreal.getTransitions()(0, 0);
  // Turbomole 7.5.1 Feb 2023
  reference = 0.1352004022564944;
  EXPECT_LE(result - reference, 1e-4);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> real_triplet({sys});
  real_triplet.settings.scfstab = Options::STABILITY_ANALYSIS::REAL;
  real_triplet.settings.densFitJ = Options::DENS_FITS::NONE;
  real_triplet.settings.triplet = true;
  real_triplet.run();
  result = real_triplet.getTransitions()(0, 0);
  // Turbomole 7.5.1 Feb 2023
  reference = 0.09958875581732934;
  EXPECT_LE(result - reference, 1e-4);

  sys->setSpin(2);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> spinflip({sys});
  spinflip.settings.scfstab = Options::STABILITY_ANALYSIS::SPINFLIP;
  spinflip.settings.densFitJ = Options::DENS_FITS::NONE;
  spinflip.run();
  result = spinflip.getTransitions()(0, 0);
  // Orca 5.0.3 Feb 2023
  reference = -0.106860;
  EXPECT_LE(result - reference, 1e-4);

  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);
  sys->setSpin(-2);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> spinflip2({sys});
  spinflip2.settings.scfstab = Options::STABILITY_ANALYSIS::SPINFLIP;
  spinflip2.settings.densFitJ = Options::DENS_FITS::NONE;
  spinflip2.run();
  result = spinflip2.getTransitions()(0, 0);
  EXPECT_LE(result - reference, 1e-4);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(LRSCFTaskTDDFTTest, TDDFT_SCFStabRootFollowing) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings settings = sys->getSettings();
  settings.basis.label = "DEF2-SVP";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE0;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  // TDDFT
  LRSCFTask<Options::SCF_MODES::RESTRICTED> real({sys});
  real.settings.scfstab = Options::STABILITY_ANALYSIS::REAL;
  real.settings.densFitJ = Options::DENS_FITS::NONE;
  real.settings.stabroot = 1;
  real.settings.stabscal = 0.8;
  real.run();
  double result = real.getTransitions()(0, 0);
  // Turbomole 7.5.1 Feb 2023
  double reference = 0.1603296108583212;
  EXPECT_LE(std::abs(result - reference), 1e-4);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
