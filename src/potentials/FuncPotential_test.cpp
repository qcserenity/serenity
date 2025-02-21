/**
 * @file FuncPotential_test.cpp
 * @author: Kevin Klahr
 *
 * @date 1. December 2016
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
#include "potentials/FuncPotential.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/FockMatrix.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "potentials/Potential.h"
#include "settings/DFTOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class FuncPotentialTest : public ::testing::Test {
 protected:
  FuncPotentialTest() {
  }

  virtual ~FuncPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_LDA) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = funcPot.getMatrix();

  Eigen::MatrixXd Fref(F.rows(), F.cols());
  Fref << -0.57441268, -0.35994693, -0.18271077, 0.00000000, -0.00000000, -0.04902895, -0.16012787, -0.22667697,
      -0.15485856, 0.00000000, -0.00000000, 0.31414675, -0.35994693, -0.39073307, -0.26090443, -0.00000000, -0.00000000,
      -0.06919500, -0.22667697, -0.30277168, -0.23298247, -0.00000000, -0.00000000, 0.21569841, -0.18271077,
      -0.26090443, -0.22635598, -0.00000000, -0.00000000, -0.04509360, -0.15485856, -0.23298247, -0.21254820,
      0.00000000, -0.00000000, 0.08411795, -0.00000000, -0.00000000, -0.00000000, -0.41596543, -0.00000000, -0.00000000,
      0.00000000, 0.00000000, 0.00000000, -0.20426443, 0.00000000, 0.00000000, -0.00000000, -0.00000000, -0.00000000,
      -0.00000000, -0.41596543, 0.00000000, 0.00000000, -0.00000000, -0.00000000, -0.00000000, -0.20426443, -0.00000000,
      -0.04902895, -0.06919500, -0.04509360, -0.00000000, 0.00000000, -0.45412511, -0.31414675, -0.21569841,
      -0.08411795, -0.00000000, -0.00000000, 0.16418006, -0.16012787, -0.22667697, -0.15485856, 0.00000000, 0.00000000,
      -0.31414675, -0.57441268, -0.35994693, -0.18271077, 0.00000000, 0.00000000, 0.04902895, -0.22667697, -0.30277168,
      -0.23298247, 0.00000000, -0.00000000, -0.21569841, -0.35994693, -0.39073307, -0.26090443, -0.00000000, 0.00000000,
      0.06919500, -0.15485856, -0.23298247, -0.21254820, 0.00000000, -0.00000000, -0.08411795, -0.18271077, -0.26090443,
      -0.22635598, 0.00000000, 0.00000000, 0.04509360, 0.00000000, -0.00000000, 0.00000000, -0.20426443, 0.00000000,
      -0.00000000, 0.00000000, -0.00000000, 0.00000000, -0.41596543, 0.00000000, 0.00000000, -0.00000000, -0.00000000,
      -0.00000000, 0.00000000, -0.20426443, -0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, -0.41596543,
      0.00000000, 0.31414675, 0.21569841, 0.08411795, 0.00000000, -0.00000000, 0.16418006, 0.04902895, 0.06919500,
      0.04509360, 0.00000000, 0.00000000, -0.45412511;

  for (unsigned i = 0; i < F.rows(); i++) {
    for (unsigned j = 0; j < F.cols(); j++) {
      EXPECT_NEAR(F(i, j), Fref(i, j), 1e-6);
    }
  }
}

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_LDA_UNRES) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(funcPot.getMatrix()));

  SPMatrix<Options::SCF_MODES::RESTRICTED> Fref(F.alpha.rows(), F.alpha.rows());

  Fref << -0.57441269, -0.35994693, -0.18271077, -0.00000000, -0.00000000, -0.04902895, -0.16012787, -0.22667697,
      -0.15485856, 0.00000000, 0.00000000, 0.31414676, -0.35994693, -0.39073307, -0.26090443, 0.00000000, 0.00000000,
      -0.06919500, -0.22667697, -0.30277168, -0.23298246, 0.00000000, 0.00000000, 0.21569841, -0.18271077, -0.26090443,
      -0.22635598, -0.00000000, -0.00000000, -0.04509360, -0.15485856, -0.23298246, -0.21254820, 0.00000000, 0.00000000,
      0.08411795, -0.00000000, -0.00000000, -0.00000000, -0.41596543, 0.00000000, -0.00000000, 0.00000000, -0.00000000,
      0.00000000, -0.20426443, -0.00000000, -0.00000000, -0.00000000, 0.00000000, 0.00000000, 0.00000000, -0.41596543,
      0.00000000, -0.00000000, 0.00000000, -0.00000000, -0.00000000, -0.20426443, -0.00000000, -0.04902895, -0.06919500,
      -0.04509360, 0.00000000, 0.00000000, -0.45412511, -0.31414676, -0.21569841, -0.08411795, 0.00000000, -0.00000000,
      0.16418006, -0.16012787, -0.22667697, -0.15485856, 0.00000000, -0.00000000, -0.31414676, -0.57441269, -0.35994693,
      -0.18271077, 0.00000000, -0.00000000, 0.04902895, -0.22667697, -0.30277168, -0.23298246, -0.00000000, 0.00000000,
      -0.21569841, -0.35994693, -0.39073307, -0.26090443, 0.00000000, 0.00000000, 0.06919500, -0.15485856, -0.23298246,
      -0.21254820, 0.00000000, -0.00000000, -0.08411795, -0.18271077, -0.26090443, -0.22635598, -0.00000000,
      -0.00000000, 0.04509360, 0.00000000, 0.00000000, 0.00000000, -0.20426443, -0.00000000, 0.00000000, 0.00000000,
      0.00000000, -0.00000000, -0.41596543, 0.00000000, 0.00000000, 0.00000000, 0.00000000, -0.00000000, -0.00000000,
      -0.20426443, -0.00000000, -0.00000000, 0.00000000, -0.00000000, 0.00000000, -0.41596543, 0.00000000, 0.31414676,
      0.21569841, 0.08411795, -0.00000000, -0.00000000, 0.16418006, 0.04902895, 0.06919500, 0.04509360, 0.00000000,
      0.00000000, -0.45412511;

  for (unsigned i = 0; i < F.rows(); i++) {
    for (unsigned j = 0; j < F.cols(); j++) {
      for_spin(F) {
        EXPECT_NEAR(F_spin(i, j), Fref(i, j), 1e-6);
      };
    }
  }
}

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_GGA) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = funcPot.getMatrix();
  Eigen::MatrixXd Fref(F.rows(), F.cols());
  Fref << -0.61072924, -0.37691231, -0.19074153, -0.00000000, -0.00000000, -0.05083544, -0.16794150, -0.23708736,
      -0.16166006, -0.00000000, -0.00000000, 0.32982024, -0.37691231, -0.40248997, -0.26694842, 0.00000000, -0.00000000,
      -0.07272697, -0.23708736, -0.31224301, -0.23861279, -0.00000000, 0.00000000, 0.22395203, -0.19074153, -0.26694842,
      -0.22783477, -0.00000000, -0.00000000, -0.04742038, -0.16166006, -0.23861279, -0.21417415, -0.00000000,
      0.00000000, 0.08747430, -0.00000000, 0.00000000, -0.00000000, -0.42730841, -0.00000000, -0.00000000, -0.00000000,
      -0.00000000, -0.00000000, -0.20982739, 0.00000000, -0.00000000, -0.00000000, -0.00000000, -0.00000000,
      -0.00000000, -0.42730841, 0.00000000, -0.00000000, -0.00000000, 0.00000000, 0.00000000, -0.20982739, 0.00000000,
      -0.05083544, -0.07272697, -0.04742038, -0.00000000, 0.00000000, -0.46998471, -0.32982024, -0.22395203,
      -0.08747430, 0.00000000, 0.00000000, 0.16971738, -0.16794150, -0.23708736, -0.16166006, -0.00000000, -0.00000000,
      -0.32982024, -0.61072924, -0.37691231, -0.19074153, -0.00000000, 0.00000000, 0.05083544, -0.23708736, -0.31224301,
      -0.23861279, -0.00000000, -0.00000000, -0.22395203, -0.37691231, -0.40248997, -0.26694842, 0.00000000, 0.00000000,
      0.07272697, -0.16166006, -0.23861279, -0.21417415, -0.00000000, 0.00000000, -0.08747430, -0.19074153, -0.26694842,
      -0.22783477, 0.00000000, -0.00000000, 0.04742038, -0.00000000, -0.00000000, -0.00000000, -0.20982739, 0.00000000,
      0.00000000, -0.00000000, 0.00000000, 0.00000000, -0.42730841, -0.00000000, 0.00000000, -0.00000000, 0.00000000,
      0.00000000, 0.00000000, -0.20982739, 0.00000000, 0.00000000, 0.00000000, -0.00000000, -0.00000000, -0.42730841,
      0.00000000, 0.32982024, 0.22395203, 0.08747430, -0.00000000, 0.00000000, 0.16971738, 0.05083544, 0.07272697,
      0.04742038, 0.00000000, 0.00000000, -0.46998471;

  for (unsigned i = 0; i < F.rows(); i++) {
    for (unsigned j = 0; j < F.cols(); j++) {
      EXPECT_NEAR(F(i, j), Fref(i, j), 1e-6);
    }
  }
}

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_GGA_UNRES) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(funcPot.getMatrix()));
  SPMatrix<RESTRICTED> Fref(F.alpha.rows(), F.alpha.rows());

  Fref << -0.61072925, -0.37691231, -0.19074153, 0.00000000, 0.00000000, -0.05083544, -0.16794150, -0.23708736, -0.16166006,
      0.00000000, -0.00000000, 0.32982025, -0.37691231, -0.40248996, -0.26694841, 0.00000000, -0.00000000, -0.07272697,
      -0.23708736, -0.31224301, -0.23861279, -0.00000000, -0.00000000, 0.22395203, -0.19074153, -0.26694841, -0.22783476,
      0.00000000, 0.00000000, -0.04742038, -0.16166006, -0.23861279, -0.21417414, 0.00000000, 0.00000000, 0.08747430,
      0.00000000, 0.00000000, 0.00000000, -0.42730840, -0.00000000, -0.00000000, 0.00000000, -0.00000000, 0.00000000,
      -0.20982739, 0.00000000, -0.00000000, 0.00000000, -0.00000000, 0.00000000, -0.00000000, -0.42730840, -0.00000000,
      -0.00000000, 0.00000000, -0.00000000, -0.00000000, -0.20982739, 0.00000000, -0.05083544, -0.07272697, -0.04742038,
      -0.00000000, -0.00000000, -0.46998471, -0.32982025, -0.22395203, -0.08747430, -0.00000000, 0.00000000, 0.16971738,
      -0.16794150, -0.23708736, -0.16166006, 0.00000000, -0.00000000, -0.32982025, -0.61072925, -0.37691231, -0.19074153,
      -0.00000000, 0.00000000, 0.05083544, -0.23708736, -0.31224301, -0.23861279, -0.00000000, 0.00000000, -0.22395203,
      -0.37691231, -0.40248996, -0.26694841, -0.00000000, 0.00000000, 0.07272697, -0.16166006, -0.23861279, -0.21417414,
      0.00000000, -0.00000000, -0.08747430, -0.19074153, -0.26694841, -0.22783476, -0.00000000, -0.00000000, 0.04742038,
      0.00000000, -0.00000000, 0.00000000, -0.20982739, -0.00000000, -0.00000000, -0.00000000, -0.00000000, -0.00000000,
      -0.42730840, 0.00000000, -0.00000000, -0.00000000, -0.00000000, 0.00000000, 0.00000000, -0.20982739, 0.00000000,
      0.00000000, 0.00000000, -0.00000000, 0.00000000, -0.42730840, 0.00000000, 0.32982025, 0.22395203, 0.08747430,
      -0.00000000, 0.00000000, 0.16971738, 0.05083544, 0.07272697, 0.04742038, -0.00000000, 0.00000000, -0.46998471;

  for (unsigned i = 0; i < F.rows(); i++) {
    for (unsigned j = 0; j < F.cols(); j++) {
      for_spin(F) {
        EXPECT_NEAR(F_spin(i, j), Fref(i, j), 1e-6);
      };
    }
  }
}

/**
 * @test FuncPotentialTest
 * @brief Tests the LDA XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_LDA) {
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));
  auto result = funcPot.getGeomGradients();

  EXPECT_NEAR(result(0, 0), 0.0, 1e-6);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-6);
  EXPECT_NEAR(result(0, 2), -2.379533e-01, 1e-6);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-6);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-6);
  EXPECT_NEAR(result(1, 2), 2.379533e-01, 1e-6);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the unrestricted LDA XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_LDA_UNRES) {
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));
  auto result = funcPot.getGeomGradients();

  EXPECT_NEAR(result(0, 0), 0.0, 1e-6);
  EXPECT_NEAR(result(0, 1), 0.0, 1e-6);
  EXPECT_NEAR(result(0, 2), -2.379533e-01, 1e-6);
  EXPECT_NEAR(result(1, 0), 0.0, 1e-6);
  EXPECT_NEAR(result(1, 1), 0.0, 1e-6);
  EXPECT_NEAR(result(1, 2), 2.379533e-01, 1e-6);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the PBE XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_GGA) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE));
  auto result = funcPot.getGeomGradients();

  EXPECT_NEAR(result(0, 0), 0.0, 5e-6);
  EXPECT_NEAR(result(0, 1), 0.0, 5e-6);
  EXPECT_NEAR(result(0, 2), -0.25023729, 5e-6);
  EXPECT_NEAR(result(1, 0), 0.0, 5e-6);
  EXPECT_NEAR(result(1, 1), 0.0, 5e-6);
  EXPECT_NEAR(result(1, 2), 0.25023729, 5e-6);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the unrestricted PBE XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_GGA_UNRES) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE_NORI);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE));
  auto result = funcPot.getGeomGradients();

  EXPECT_NEAR(result(0, 0), 0.0, 5e-6);
  EXPECT_NEAR(result(0, 1), 0.0, 5e-6);
  EXPECT_NEAR(result(0, 2), -0.25023729, 5e-6);
  EXPECT_NEAR(result(1, 0), 0.0, 5e-6);
  EXPECT_NEAR(result(1, 1), 0.0, 5e-6);
  EXPECT_NEAR(result(1, 2), 0.25023729, 5e-6);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the BP86 matrix-potential calculated 2 ways.
 */
TEST_F(FuncPotentialTest, H2_Potential_GGA) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), dMat->getDensityMatrix().getBasisController(), grid);
  basisFunctionOnGridController->setHighestDerivative(2);
  auto densOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(basisFunctionOnGridController, 0.0);
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>(densOnGridCalculator, dMat, 2);
  auto gridToMatrix =
      std::make_shared<ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED>>(basisFunctionOnGridController, 0.0);

  FunctionalLibrary<Options::SCF_MODES::RESTRICTED> flib(128);
  auto funcData1 = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS,
                                 CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86),
                                 densOnGridController);

  auto funcData2 = flib.calcData(FUNCTIONAL_DATA_TYPE::POTENTIAL,
                                 CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86),
                                 densOnGridController);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F1(dMat->getDensityMatrix().getBasisController());
  FockMatrix<Options::SCF_MODES::RESTRICTED> F2(dMat->getDensityMatrix().getBasisController());

  gridToMatrix->addScalarOperatorToMatrix(F2, *funcData2.potential);
  gridToMatrix->addScalarOperatorToMatrix(F1, *funcData1.dFdRho, *funcData1.dFdGradRho);

  for (unsigned i = 0; i < F1.size(); i++) {
    EXPECT_NEAR(F1.data()[i], F2.data()[i], 1e-5);
  }
}

TEST_F(FuncPotentialTest, CompositeVsCustomFunctional) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86, true);
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemController, dMat, grid, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.basis.label = "6-31GS";
  settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::X_B88,
                                          BasicFunctionals::BASIC_FUNCTIONALS::C_P86};
  settings.customFunc.mixingFactors = {1.0, 1.0};
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::EITHER_OR;
  auto systemController2 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86, settings, 0, 0);
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat2 =
      systemController2->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  std::shared_ptr<GridController> grid2 = systemController2->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot2(systemController2, dMat2, grid2, Functional(settings.customFunc));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = funcPot.getMatrix();
  FockMatrix<Options::SCF_MODES::RESTRICTED> F2 = funcPot2.getMatrix();

  for (unsigned i = 0; i < F.size(); i++) {
    EXPECT_NEAR(F.data()[i], F2.data()[i], 1e-8);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController2);
}

} // namespace Serenity
