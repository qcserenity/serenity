/**
 * @file   GradientTask_test.cpp
 *
 * @date   June 28, 2017
 * @author Kevin Klahr
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
#include "tasks/GradientTask.h"
#include "data/ElectronicStructure.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "math/Matrix.h"
#include "settings/Options.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>

namespace Serenity {
/**
 * @class GradientTaskTest
 * @brief Sets everything up for the tests of GradientTask.h/.cpp .
 *
 */
class GradientTaskTest : public ::testing::Test {
 protected:
  GradientTaskTest() = default;

  virtual ~GradientTaskTest() = default;

  static void SetUpTestCase() {
    auto& libint = Libint::getInstance();
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);
  }

  virtual void TearDown() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
    auto& libint = Libint::getInstance();
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 4);
  }
};

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant gradients for H2O with BP86
 */
TEST_F(GradientTaskTest, D3_Grad_BP86) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  Settings settings = _activeSystem->getSettings();
  settings.dft.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
  _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem});
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  // 14.09.2020 Changed due to switch from full four center integral based gradients to RI.
  EXPECT_NEAR(actGradient(0, 0), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.0068987965582244081, 4e-5);
  EXPECT_NEAR(actGradient(1, 0), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), -0.0068987965582246302, 4e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant FDE gradients for the H2O dimer with BP86
 */
TEST_F(GradientTaskTest, D3_Grad_FDE_BP86) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem}, {_environmentSystem});
  task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.049598237810465171, 4e-5);
  EXPECT_NEAR(actGradient(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), 0.55318994972444768, 4e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant gradients for H2O with BP86
 */
TEST_F(GradientTaskTest, TranslationalInvariance_BP86) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem});
  task.settings.transInvar = true;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0, 0), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), 0.0034355317315704157, 2e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.0034355317315703914, 2e-5);
  EXPECT_NEAR(actGradient(1, 0), -0.0034355317315703914, 2e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), -0.0034355317315704153, 2e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant FDE gradients for the H2O dimer with BP86
 */
TEST_F(GradientTaskTest, TranslationalInvarianceFDE_BP86) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem}, {_environmentSystem});
  task.settings.transInvar = true;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), -0.27661202908904003, 2e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.024780391604813233, 2e-5);
  EXPECT_NEAR(actGradient(1, 0), -0.02478039160481324, 2e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), 0.27661202908903998, 2e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of a Ne dimer in a 6-31Gs basis.
 */
TEST_F(GradientTaskTest, NeDimerGradient) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;

  // Neon dimer
  auto systemControllerNe = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs);
  auto es = systemControllerNe->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials;
  potentials =
      systemControllerNe->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerNe});
  task.settings.transInvar = false;
  task.settings.print = false;
  task.run();

  Eigen::MatrixXd gradientsNe = systemControllerNe->getGeometry()->getGradients();

  EXPECT_NEAR(gradientsNe(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(gradientsNe(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(gradientsNe(0, 2), 0.033594406857773151, 1e-5);
  EXPECT_NEAR(gradientsNe(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(gradientsNe(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(gradientsNe(1, 2), -0.03359440685778825, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of F2 in a 6-31Gs basis.
 */
TEST_F(GradientTaskTest, F2Gradient) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;

  // F2
  auto systemControllerF2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs, true);
  auto es = systemControllerF2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials;
  potentials =
      systemControllerF2->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerF2});
  task.settings.transInvar = false;
  task.settings.print = false;
  task.run();

  Eigen::MatrixXd gradientsF2 = systemControllerF2->getGeometry()->getGradients();

  EXPECT_NEAR(gradientsF2(0, 0), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(0, 1), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(0, 2), -0.1259425, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 0), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 1), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 2), 0.1259425, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of F2 in a 6-31Gs basis numerically.
 */

TEST_F(GradientTaskTest, F2GradientsNum) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;

  // F2
  auto systemControllerF2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs, true);
  auto es = systemControllerF2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials;
  potentials =
      systemControllerF2->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerF2});
  task.settings.transInvar = false;
  task.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
  task.settings.print = false;
  task.run();

  Eigen::MatrixXd gradientsF2 = systemControllerF2->getGeometry()->getGradients();

  EXPECT_NEAR(gradientsF2(0, 0), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(0, 1), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(0, 2), -0.1259425, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 0), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 1), 0.0000000, 1e-5);
  EXPECT_NEAR(gradientsF2(1, 2), 0.1259425, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of F2 in a 6-31Gs basis + PCM numerical vs analytical.
 */

TEST_F(GradientTaskTest, F2CPCMGradients) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.pcm.use = true;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::CPCM;
  settings.scf.diisThreshold = 1e-7;
  settings.scf.energyThreshold = 1e-8;
  settings.scf.rmsdThreshold = 1e-8;

  // F2
  auto systemControllerF2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs, settings);
  auto es = systemControllerF2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials;
  potentials =
      systemControllerF2->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> numTask({systemControllerF2});
  numTask.settings.transInvar = false;
  numTask.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
  numTask.settings.print = false;
  numTask.settings.numGradStepSize = 0.01;
  numTask.run();

  Eigen::MatrixXd numGrads = systemControllerF2->getGeometry()->getGradients();

  GradientTask<Options::SCF_MODES::RESTRICTED> analyTask({systemControllerF2});
  analyTask.settings.transInvar = false;
  analyTask.settings.gradType = Options::GRADIENT_TYPES::ANALYTICAL;
  analyTask.run();

  Eigen::MatrixXd analyGrads = systemControllerF2->getGeometry()->getGradients();
  // Note that analytical and numerical CPCM gradients deviate somewhat, since the
  // change in tesserae area is not respected in the analytical gradients.
  EXPECT_NEAR(analyGrads(0, 0), numGrads(0, 0), 5e-5);
  EXPECT_NEAR(analyGrads(0, 1), numGrads(0, 1), 5e-5);
  EXPECT_NEAR(analyGrads(0, 2), numGrads(0, 2), 5e-5);
  EXPECT_NEAR(analyGrads(1, 0), numGrads(1, 0), 5e-5);
  EXPECT_NEAR(analyGrads(1, 1), numGrads(1, 1), 5e-5);
  EXPECT_NEAR(analyGrads(1, 2), numGrads(1, 2), 5e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemControllerF2);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test GradientTaskTest
 * @brief Tests the FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDA) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController>> allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials1;
  potentials1 = _activeSystem->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials2;
  potentials2 =
      _environmentSystem->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem}, {_environmentSystem});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.049560783209758874, 1e-3);
  EXPECT_NEAR(actGradient(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), 0.55322405817816644, 1e-3);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the numerical FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDANum) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE, true);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController>> allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials1;
  potentials1 = _activeSystem->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> potentials2;
  potentials2 =
      _environmentSystem->getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  {
    GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem}, {_environmentSystem});
    task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
    task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
    task.settings.FDEgridCutOff = -1.0;
    task.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
    task.settings.print = false;
    task.run();
  }

  auto actGradientNum = _activeSystem->getGeometry()->getGradients();

  {
    GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem}, {_environmentSystem});
    task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
    task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
    task.settings.FDEgridCutOff = -1.0;
    task.settings.gradType = Options::GRADIENT_TYPES::ANALYTICAL;
    task.settings.print = false;
    task.run();
  }

  auto actGradientAna = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradientNum(0, 0), actGradientAna(0, 0), 1e-5);
  EXPECT_NEAR(actGradientNum(0, 1), actGradientAna(0, 1), 1e-5);
  EXPECT_NEAR(actGradientNum(0, 2), actGradientAna(0, 2), 1e-3);
  EXPECT_NEAR(actGradientNum(1, 0), actGradientAna(1, 0), 1e-5);
  EXPECT_NEAR(actGradientNum(1, 1), actGradientAna(1, 1), 1e-5);
  EXPECT_NEAR(actGradientNum(1, 2), actGradientAna(1, 2), 1e-3);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDA_UNRES) {
  auto _activeSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController>> allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED>> potentials1;
  potentials1 =
      _activeSystem->getPotentials<Options::SCF_MODES::UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED>> potentials2;
  potentials2 =
      _environmentSystem->getPotentials<Options::SCF_MODES::UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  GradientTask<Options::SCF_MODES::UNRESTRICTED> task({_activeSystem}, {_environmentSystem});
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::TF;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0, 2), 0.041510210877332976, 1e-3);
  EXPECT_NEAR(actGradient(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1, 2), 0.60553010482453451, 1e-3);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

TEST_F(GradientTaskTest, FDEGradientsVsNumerical) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Def2_TZVP_DFT, true);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He_Def2_TZVP_DFT, true);

  GradientTask<Options::SCF_MODES::RESTRICTED> analytical({a, b}, {});
  analytical.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  analytical.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  analytical.settings.FDEgridCutOff = -1.0;
  analytical.settings.print = false;
  analytical.settings.gradType = Options::GRADIENT_TYPES::ANALYTICAL;
  analytical.run();

  Matrix<double> aGradAnalytical = a->getGeometry()->getGradients();
  Matrix<double> bGradAnalytical = b->getGeometry()->getGradients();

  GradientTask<Options::SCF_MODES::RESTRICTED> numerical({a, b}, {});
  numerical.settings = analytical.settings;
  numerical.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
  numerical.settings.numGradStepSize = 0.001;
  numerical.run();

  Matrix<double> aGradNumerical = a->getGeometry()->getGradients();
  Matrix<double> bGradNumerical = b->getGeometry()->getGradients();

  double aDiff = (aGradNumerical - aGradAnalytical).array().abs().maxCoeff();
  double bDiff = (bGradNumerical - bGradAnalytical).array().abs().maxCoeff();
  // The gradient should be in the order of magnitude of 1e-2 to 1e-3. We should at least get the first few digits right.
  EXPECT_NEAR(aDiff, 0.0, 1e-4);
  EXPECT_NEAR(bDiff, 0.0, 1e-4);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

TEST_F(GradientTaskTest, ElectricField) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, true);
  Settings settings = _activeSystem->getSettings();
  settings.efield.use = true;
  settings.efield.analytical = true;
  settings.efield.pos2 = {1, 1, 1};
  settings.efield.fieldStrength = 1e-3;
  _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, settings);

  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem});
  // task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  // TURBOMOLE 7.5.1 Oct 2021.
  Eigen::MatrixXd reference = actGradient;
  reference.row(0) << 0.15059558720612e-02, -.65188084013720e-01, 0.52813351777630e-02;
  reference.row(1) << -.31067293123075e-01, 0.43324559589465e-02, -.19340843317178e-02;
  reference.row(2) << 0.24041258523979e-01, 0.61156419853363e-01, -.34626500435051e-02;
  reference.row(3) << 0.16194159680843e-01, 0.14551320831473e-02, -.71917926741680e-01;
  reference.row(4) << 0.10662968508014e-01, 0.17457747113100e-01, 0.15091892929620e-01;
  reference.row(5) << -.21337049462470e-01, -.19213670994811e-01, 0.56941433009442e-01;

  double maxDiff = (reference - actGradient).cwiseAbs().maxCoeff();
  double accuracy = 1e-7;

  EXPECT_LE(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(_activeSystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(GradientTaskTest, PointChargeGradients) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::ETHANOL_def2_SVP_HF);
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/" + sys->getSystemName() + "/point-charges.pc";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  Settings settings = sys->getSettings();
  settings.extCharges.externalChargesFile = pathToTestsResources;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::ETHANOL_def2_SVP_HF, settings);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.run();

  const auto& gradients = sys->getGeometry()->getGradients();
  const auto& pointChargeGradients = sys->getPointChargeGradients();

  // Reference generated with Orca 4.2.0
  Eigen::MatrixXd referenceGradients(9, 3);
  referenceGradients << 0.00918164, -0.00930644, 0.00432299, -0.0104084, 0.0270658, -0.00271828, -0.00519274,
      -0.00452777, 0.00452071, -0.00225914, -0.00382822, -0.000416061, -0.00231768, 0.00308402, -0.00529075, 0.00163825,
      -0.00353335, 0.0278773, 0.00215748, -0.000452374, -0.00555908, -0.00920128, -0.00455667, 0.00422912, 0.0164212,
      -0.00409334, -0.027123;
  Eigen::MatrixXd referencePointChargeGradients(10, 3);
  referencePointChargeGradients << 1.54667e-06, 2.28381e-06, 3.23415e-06, -7.55103e-07, -9.68785e-07, -1.51259e-06,
      -8.73348e-07, -1.13626e-06, -1.76019e-06, -6.80872e-06, 1.2159e-08, -6.83198e-07, 3.37468e-06, 1.05529e-07,
      2.11097e-07, 3.10415e-06, -3.4387e-08, 3.34131e-07, -2.39119e-06, 3.21388e-06, -2.39693e-07, 1.292e-06,
      -1.59313e-06, 1.38987e-07, 1.09579e-06, -1.47915e-06, 1.03665e-07, -9.37566e-06, 2.414e-05, -1.02289e-05;

  const double diffGrad = (gradients - referenceGradients).array().abs().maxCoeff();
  // The reference gradients are only provided up to the 6th digit.
  EXPECT_NEAR(diffGrad, 0.0, 1e-5);
  const double pointChargeDiff = (pointChargeGradients.topRows(10) - referencePointChargeGradients).array().abs().maxCoeff();
  // The point charge gradients are smaller and available up to a higher accuracy.
  EXPECT_NEAR(pointChargeDiff, 0.0, 1e-8);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
