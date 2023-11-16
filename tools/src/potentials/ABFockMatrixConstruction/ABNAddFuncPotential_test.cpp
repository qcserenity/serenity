/**
 * @file ABNAddFuncPotential_test.cpp
 *
 * @date May 18, 2018
 * @author: Moritz Bensberg
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
#include "potentials/ABFockMatrixConstruction/ABNAddFuncPotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "potentials/NAddFuncPotential.h"
#include "settings/DFTOptions.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class ABNAddFuncPotentialTest ABCoulombInteractionPotential_test.cpp
 * @brief Sets up everything for the test of ABNAddFuncPotential.h/.cpp.
 */
class ABNAddFuncPotentialTest : public ::testing::Test {
 protected:
  ABNAddFuncPotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~ABNAddFuncPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
};

/**
 * @test ABNAddFuncPotentialTest
 * @brief Tests the LDA AA Fock Matrix of an H2 dimer. Identical to the equivalent
 *        test in NAddFuncPotentialTest.
 */
TEST_F(ABNAddFuncPotentialTest, H2_dimer_rFockMatrixAAB_LDA) {
  auto basisA = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerB->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABNAddFuncPotential<Options::SCF_MODES::RESTRICTED> abNAddPot(
      systemControllerA, basisA, basisA, envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  NAddFuncPotential<Options::SCF_MODES::RESTRICTED> nAddPot(
      systemControllerA,
      systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(),
      envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F_AB = abNAddPot.getMatrix();
  SPMatrix<Options::SCF_MODES::RESTRICTED> F_AA = nAddPot.getMatrix();

  auto diff = (F_AB - F_AA).array().abs().maxCoeff();
  // accept only white noise!
  EXPECT_NEAR(diff, 0.0, 1e-12);
}

/**
 * @test ABNAddFuncPotentialTest
 * @brief Tests the LDA AB Fock Matrix of an H2 dimer.
 */
TEST_F(ABNAddFuncPotentialTest, H2_dimer_rFockMatrixABB_LDA) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerB->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABNAddFuncPotential<Options::SCF_MODES::RESTRICTED> abNAddPot(
      systemControllerA, basisA, basisB, envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F = abNAddPot.getMatrix();

  EXPECT_NEAR(F(0, 0), -0.0083840258804812658, 1e-5);
  EXPECT_NEAR(F(1, 0), -0.055300758004326538, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.12831003302437624, 1e-5);
  EXPECT_NEAR(F(0, 1), -0.005639132369059649, 1e-5);
  EXPECT_NEAR(F(0, 2), -0.0039221790739654878, 1e-5);
  EXPECT_NEAR(F(0, 3), -0.0044904816299493342, 1e-5);
}
/**
 * @test ABNAddFuncPotentialTest
 * @brief Tests the GGA AB Fock Matrix of an H2 dimer.
 */
TEST_F(ABNAddFuncPotentialTest, H2_dimer_rFockMatrixABB_GGA) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerB->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABNAddFuncPotential<Options::SCF_MODES::RESTRICTED> abNAddPot(
      systemControllerA, basisA, basisB, envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F = abNAddPot.getMatrix();

  EXPECT_NEAR(F(0, 0), -0.0076542762861270645, 1e-5);
  EXPECT_NEAR(F(1, 0), -0.055253744527196447, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.12698015387994177, 1e-5);
  EXPECT_NEAR(F(0, 1), -0.0047051305109973322, 1e-5);
  EXPECT_NEAR(F(0, 2), -0.0037534727526342374, 1e-5);
  EXPECT_NEAR(F(0, 3), -0.0038349418254594093, 1e-5);
}
/**
 * @test ABNAddFuncPotentialTest
 * @brief Tests the GGA AB Fock Matrix of an H2 dimer for multiple environment
 *        density matrices.
 */
TEST_F(ABNAddFuncPotentialTest, H2_dimer_rFockMatrixABBB_GGA) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      2, systemControllerB->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABNAddFuncPotential<Options::SCF_MODES::RESTRICTED> abNAddPot(
      systemControllerA, basisA, basisB, envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F = abNAddPot.getMatrix();

  EXPECT_NEAR(F(0, 0), -0.012352511055803322, 1e-5);
  EXPECT_NEAR(F(1, 0), -0.080462987092342453, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.18794150046341079, 1e-5);
  EXPECT_NEAR(F(0, 1), -0.0081584035257395351, 1e-5);
  EXPECT_NEAR(F(0, 2), -0.0057641379450071991, 1e-5);
  EXPECT_NEAR(F(0, 3), -0.0065285874827300102, 1e-5);
}
/**
 * @test ABNAddFuncPotentialTest
 * @brief Tests the unrestricted GGA AB Fock Matrix of an H2 dimer.
 */
TEST_F(ABNAddFuncPotentialTest, H2_dimer_uFockMatrixABB_GGA) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> envDmat(
      1, systemControllerB->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());

  ABNAddFuncPotential<Options::SCF_MODES::UNRESTRICTED> abNAddPot(
      systemControllerA, basisA, basisB, envDmat, systemControllerA->getGridController(),
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  SPMatrix<Options::SCF_MODES::UNRESTRICTED> F = abNAddPot.getMatrix();

  EXPECT_NEAR(F.alpha(0, 0), -0.0068480563425059567, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), -0.050840539961720332, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), -0.11455604883352168, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), -0.0044249629657753245, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), -0.0033853933966209167, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), -0.0035719792926375772, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), -0.0068480563425059567, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), -0.050840539961720332, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), -0.11455604883352168, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), -0.0044249629657753245, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), -0.0033853933966209167, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), -0.0035719792926375772, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

} /* namespace Serenity */
