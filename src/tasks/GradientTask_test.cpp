/**
 * @file   GradientTask_test.cpp
 *
 * @date   June 28, 2017
 * @author Kevin Klahr
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
#include "data/ElectronicStructure.h"
#include "tasks/GradientTask.h"
#include "settings/Options.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <cmath>
#include <gtest/gtest.h>

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
    libint.keepEngines(libint2::Operator::coulomb,0,2);
    libint.keepEngines(libint2::Operator::coulomb,0,3);
    libint.keepEngines(libint2::Operator::coulomb,0,4);
    libint.keepEngines(libint2::Operator::coulomb,1,2);
    libint.keepEngines(libint2::Operator::coulomb,1,3);
    libint.keepEngines(libint2::Operator::coulomb,1,4);
  }

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
    auto& libint = Libint::getInstance();
    libint.freeEngines(libint2::Operator::coulomb,0,2);
    libint.freeEngines(libint2::Operator::coulomb,0,3);
    libint.freeEngines(libint2::Operator::coulomb,0,4);
    libint.freeEngines(libint2::Operator::coulomb,1,2);
    libint.freeEngines(libint2::Operator::coulomb,1,3);
    libint.freeEngines(libint2::Operator::coulomb,1,4);
  }

};


/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant gradients for H2O with BP86
 */
TEST_F(GradientTaskTest, D3_Grad_BP86) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto actSettings = _activeSystem->getSettings();
  actSettings.dft.functional = Options::XCFUNCTIONALS::BP86;
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem});
  task.settings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
//  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0), 0.000,  1e-5);
  EXPECT_NEAR(actGradient(0,1), 0.000,  1e-5);
  EXPECT_NEAR(actGradient(0,2), 0.006959380382, 1e-5);
  EXPECT_NEAR(actGradient(1,0), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1,1), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1,2),-0.006959380382, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant FDE gradients for the H2O dimer with BP86
 */
TEST_F(GradientTaskTest, D3_Grad_FDE_BP86) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.dispersion = Options::DFT_DISPERSION_CORRECTIONS::D3BJ;
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,1),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,2),  0.049598237810465171, 1e-5);
  EXPECT_NEAR(actGradient(1,0),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,1),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,2),  0.55318994972444768, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);


}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant gradients for H2O with BP86
 */
TEST_F(GradientTaskTest, TranslationalInvariance_BP86) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem});
  task.settings.transInvar = true;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0), 0.000,  1e-5);
  EXPECT_NEAR(actGradient(0,1), 0.0034658236626372288,  1e-5);
  EXPECT_NEAR(actGradient(0,2), 0.0034658236625880194, 1e-5);
  EXPECT_NEAR(actGradient(1,0),-0.0034658236625880199, 1e-5);
  EXPECT_NEAR(actGradient(1,1), 0.000, 1e-5);
  EXPECT_NEAR(actGradient(1,2),-0.0034658236626372292, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the translationally invariant FDE gradients for the H2O dimer with BP86
 */
TEST_F(GradientTaskTest, TranslationalInvarianceFDE_BP86) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.transInvar = true;
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,1), -0.27661202908904003, 1e-5);
  EXPECT_NEAR(actGradient(0,2),  0.024780391604813233, 1e-5);
  EXPECT_NEAR(actGradient(1,0), -0.02478039160481324, 1e-5);
  EXPECT_NEAR(actGradient(1,1),  0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,2),  0.27661202908903998, 1e-5);

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

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials;
  potentials = systemControllerNe->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerNe});
  task.settings.transInvar = false;
  task.settings.print = false;
  task.run();

  Matrix<double> gradientsNe = systemControllerNe->getGeometry()->getGradients();

  EXPECT_NEAR(gradientsNe(0,0),0.0,1e-5);
  EXPECT_NEAR(gradientsNe(0,1),0.0,1e-5);
  EXPECT_NEAR(gradientsNe(0,2), 0.033594406857773151,1e-5);
  EXPECT_NEAR(gradientsNe(1,0),0.0,1e-5);
  EXPECT_NEAR(gradientsNe(1,1),0.0,1e-5);
  EXPECT_NEAR(gradientsNe(1,2),-0.03359440685778825,1e-5);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of F2 in a 6-31Gs basis.
 */

TEST_F(GradientTaskTest, F2Gradient) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;

  // F2
  auto systemControllerF2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs);
  auto es = systemControllerF2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials;
  potentials = systemControllerF2->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerF2});
  task.settings.transInvar = false;
  task.settings.print = false;
  task.run();

  Matrix<double> gradientsF2 = systemControllerF2->getGeometry()->getGradients();


  EXPECT_NEAR(gradientsF2(0,0),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(0,1),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(0,2),-0.1252599,1e-5);
  EXPECT_NEAR(gradientsF2(1,0),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(1,1),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(1,2),0.1252599,1e-5);
}

/**
 * @test GradientTaskTest
 * @brief Tests the total gradient of F2 in a 6-31Gs basis numerically.
 */

TEST_F(GradientTaskTest, F2GradientsNum) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  // F2
  auto systemControllerF2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F2_6_31Gs);
  auto es = systemControllerF2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials;
  potentials = systemControllerF2->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  es->attachPotentials(potentials);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({systemControllerF2});
  task.settings.transInvar = false;
  task.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
  task.settings.print = false;
  task.run();

  Matrix<double> gradientsF2 = systemControllerF2->getGeometry()->getGradients();

  EXPECT_NEAR(gradientsF2(0,0),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(0,1),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(0,2),-0.1252599,1e-4);
  EXPECT_NEAR(gradientsF2(1,0),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(1,1),0.0,1e-5);
  EXPECT_NEAR(gradientsF2(1,2),0.1252599,1e-4);
}

/**
 * @test GradientTaskTest
 * @brief Tests the FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDA) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController> > allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials1;
  potentials1 = _activeSystem->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials2;
  potentials2 = _environmentSystem->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.print = false;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,2), 0.049560783209758874, 1e-3);
  EXPECT_NEAR(actGradient(1,0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,2), 0.55322405817816644, 1e-3);


  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the numerical FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDANum) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE,true);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE,true);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController> > allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials1;
  potentials1 = _activeSystem->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED> > potentials2;
  potentials2 = _environmentSystem->getPotentials<Options::SCF_MODES::RESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  {GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.gradType = Options::GRADIENT_TYPES::NUMERICAL;
  task.settings.print = false;
  task.run();}

  auto actGradientNum = _activeSystem->getGeometry()->getGradients();

  {GradientTask<Options::SCF_MODES::RESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.settings.gradType = Options::GRADIENT_TYPES::ANALYTICAL;
  task.settings.print = false;
  task.run();}

  auto actGradientAna = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradientNum(0,0), actGradientAna(0,0), 1e-5);
  EXPECT_NEAR(actGradientNum(0,1), actGradientAna(0,1), 1e-5);
  EXPECT_NEAR(actGradientNum(0,2), actGradientAna(0,2), 1e-3);
  EXPECT_NEAR(actGradientNum(1,0), actGradientAna(1,0), 1e-5);
  EXPECT_NEAR(actGradientNum(1,1), actGradientAna(1,1), 1e-5);
  EXPECT_NEAR(actGradientNum(1,2), actGradientAna(1,2), 1e-3);


  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test GradientTaskTest
 * @brief Tests the FDE gradients with BP86.
 */
TEST_F(GradientTaskTest, FDEGradientsLDA_UNRES) {
  auto _activeSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto _environmentSystem =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  /*
   * Initialize data objects and gather data
   */
  // list of all systems for easy loops
  std::vector<std::shared_ptr<SystemController> > allsystems;
  allsystems.push_back(_environmentSystem);
  allsystems.push_back(_activeSystem);

  auto es1 = _activeSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto es2 = _environmentSystem->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();

  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> > potentials1;
  potentials1 = _activeSystem->getPotentials<Options::SCF_MODES::UNRESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es1->attachPotentials(potentials1);
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::UNRESTRICTED> > potentials2;
  potentials2 = _environmentSystem->getPotentials<Options::SCF_MODES::UNRESTRICTED,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  es2->attachPotentials(potentials2);

  GradientTask<Options::SCF_MODES::UNRESTRICTED> task({_activeSystem},{_environmentSystem});
  task.settings.naddKinFunc = Options::KINFUNCTIONALS::TF;
  task.settings.naddXCFunc = Options::XCFUNCTIONALS::BP86;
  task.settings.FDEgridCutOff = -1.0;
  task.run();

  auto actGradient = _activeSystem->getGeometry()->getGradients();

  EXPECT_NEAR(actGradient(0,0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(0,2), 0.041510210877332976, 1e-3);
  EXPECT_NEAR(actGradient(1,0), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,1), 0.0, 1e-5);
  EXPECT_NEAR(actGradient(1,2), 0.60553010482453451, 1e-3);


  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}
}
