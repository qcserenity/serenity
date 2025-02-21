/**
 * @file RICoulombPotential_test.cpp
 *
 * @date Dec 1, 2016
 * @author: Kevin Klahr
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
#include "potentials/RICoulombPotential.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "potentials/Potential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class RICoulombPotentialTest : public ::testing::Test {
 protected:
  RICoulombPotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~RICoulombPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
};

/**
 * @test RICoulombPotentialTest
 * @brief Tests the Fock Matrix of an H2
 */
TEST_F(RICoulombPotentialTest, H2_rFockMatrix) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  auto basisController = systemController->getAtomCenteredBasisController();
  auto auxBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  auto& factory = RI_J_IntegralControllerFactory::getInstance();
  auto ri_j_IntController = factory.produce(basisController, auxBasisController);

  RICoulombPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController, dMat, ri_j_IntController,
                                                             basisController->getPrescreeningThreshold(), 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.58366901, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test RICoulombPotentialTest
 * @brief Tests the Fock Matrix of an H2
 */
TEST_F(RICoulombPotentialTest, H2_uFockMatrix) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  auto basisController = systemController->getAtomCenteredBasisController();
  auto auxBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  auto& factory = RI_J_IntegralControllerFactory::getInstance();
  auto ri_j_IntController = factory.produce(basisController, auxBasisController);

  RICoulombPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemController, dMat, ri_j_IntController,
                                                               basisController->getPrescreeningThreshold(), 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(coulPot.getMatrix()));

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.58366901, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), 1.58366901, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test RICoulombPotentialTest
 * @brief Tests the gradient of an H2
 */
TEST_F(RICoulombPotentialTest, H2_rGradients) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  auto basisController = systemController->getAtomCenteredBasisController();
  auto auxBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  auto& factory = RI_J_IntegralControllerFactory::getInstance();
  auto ri_j_IntController = factory.produce(basisController, auxBasisController);

  RICoulombPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController, dMat, ri_j_IntController,
                                                             basisController->getPrescreeningThreshold(), 0.0, 0.0, 0);

  auto ElEl = coulPot.getGeomGradients();

  EXPECT_NEAR(ElEl(0, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(0, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(0, 2), 0.69353092420185058, 1e-4);
  EXPECT_NEAR(ElEl(1, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(1, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(1, 2), -0.69353092420185058, 1e-4);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test RICoulombPotentialTest
 * @brief Tests the gradient of an H2
 */
TEST_F(RICoulombPotentialTest, H2_uGradients) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  auto basisController = systemController->getAtomCenteredBasisController();
  auto auxBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  auto& factory = RI_J_IntegralControllerFactory::getInstance();
  auto ri_j_IntController = factory.produce(basisController, auxBasisController);

  RICoulombPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemController, dMat, ri_j_IntController,
                                                               basisController->getPrescreeningThreshold(), 0.0, 0.0, 0);

  auto ElEl = coulPot.getGeomGradients();

  EXPECT_NEAR(ElEl(0, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(0, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(0, 2), 0.69353092420185058, 1e-4);
  EXPECT_NEAR(ElEl(1, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(1, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElEl(1, 2), -0.69353092420185058, 1e-4);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

} // namespace Serenity
