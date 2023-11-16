/**
 * @file ABExchangePotential_test.cpp
 *
 * @date Jun 21, 2018
 * @author Moritz Bensberg
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
#include "potentials/ABFockMatrixConstruction/ABExchangePotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "potentials/ERIPotential.h"
#include "potentials/ExchangeInteractionPotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABExchangePotentialTest : public ::testing::Test {
 protected:
  ABExchangePotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)),
      systemControllerC(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~ABExchangePotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
  std::shared_ptr<SystemController> systemControllerC;
};

/**
 * @test
 * @brief Tests the restricted construction of the f_AA matrix for the density of A.
 */
TEST_F(ABExchangePotentialTest, H2_dimer_rFockMatrixAAA) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto auxBasisA = systemControllerA->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  //  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerA->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABExchangePotential<SPIN> abExchangePot(systemControllerA, basisA, basisA, envDmat, 1.0);

  ABCoulombInteractionPotential<SPIN> abCoulombIntPot(systemControllerA, basisA, basisA, envDmat, false,
                                                      Options::DENS_FITS::RI, auxBasisA);

  ERIPotential<SPIN> hfPotential(systemControllerA, envDmat[0], 1.0, 1E-10, 0.0, 0.0, 0, false);

  SPMatrix<SPIN> f_AAA1 = abExchangePot.getMatrix() + abCoulombIntPot.getMatrix();
  SPMatrix<SPIN> f_AAA2 = hfPotential.getMatrix();

  SPMatrix<SPIN> test = f_AAA1.transpose();

  EXPECT_NEAR((f_AAA1 - test).array().abs().maxCoeff(), 0.0, 1e-12);
  EXPECT_NEAR(f_AAA1(0, 0), f_AAA2(0, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(1, 0), f_AAA2(1, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(2, 0), f_AAA2(2, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 1), f_AAA2(0, 1), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 2), f_AAA2(0, 2), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 3), f_AAA2(0, 3), 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AA matrix for the density of C.
 */
TEST_F(ABExchangePotentialTest, H2_dimer_rFockMatrixAAC) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  //  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerC->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABExchangePotential<SPIN> abExchangePot(systemControllerA, basisA, basisA, envDmat, 1.0);

  ExchangeInteractionPotential<SPIN> exchangeInteractionPotential(basisA, envDmat, 1.0, 1E-10);

  SPMatrix<SPIN> f_AAA1 = abExchangePot.getMatrix();
  SPMatrix<SPIN> f_AAA2 = exchangeInteractionPotential.getMatrix();

  SPMatrix<SPIN> test = f_AAA1.transpose();

  EXPECT_NEAR((f_AAA1 - test).sum(), 0.0, 1e-10);
  EXPECT_NEAR(f_AAA1(0, 0), f_AAA2(0, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(1, 0), f_AAA2(1, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(2, 0), f_AAA2(2, 0), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 1), f_AAA2(0, 1), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 2), f_AAA2(0, 2), 1e-5);
  EXPECT_NEAR(f_AAA1(0, 3), f_AAA2(0, 3), 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AA matrix for the density of C.
 */
TEST_F(ABExchangePotentialTest, H2_dimer_uFockMatrixAAC) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  //  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerC->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABExchangePotential<SPIN> abExchangePot(systemControllerA, basisA, basisA, envDmat, 1.0);

  ExchangeInteractionPotential<SPIN> exchangeInteractionPotential(basisA, envDmat, 1.0, 1E-10);

  SPMatrix<SPIN> f_AAA1 = abExchangePot.getMatrix();
  SPMatrix<SPIN> f_AAA2 = exchangeInteractionPotential.getMatrix();

  double test = 0.0;
  for_spin(f_AAA1) {
    test += (f_AAA1_spin - f_AAA1_spin.transpose()).sum();
  };

  EXPECT_NEAR(test, 0.0, 1e-10);

  EXPECT_NEAR(f_AAA1.alpha(0, 0), f_AAA2.alpha(0, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(1, 0), f_AAA2.alpha(1, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(2, 0), f_AAA2.alpha(2, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 1), f_AAA2.alpha(0, 1), 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 2), f_AAA2.alpha(0, 2), 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 3), f_AAA2.alpha(0, 3), 1e-5);

  EXPECT_NEAR(f_AAA1.beta(0, 0), f_AAA2.beta(0, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.beta(1, 0), f_AAA2.beta(1, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.beta(2, 0), f_AAA2.beta(2, 0), 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 1), f_AAA2.beta(0, 1), 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 2), f_AAA2.beta(0, 2), 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 3), f_AAA2.beta(0, 3), 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AB matrix for the density of C.
 */
TEST_F(ABExchangePotentialTest, H2_dimer_rFockMatrixABC) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerC->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABExchangePotential<SPIN> abExchangePot(systemControllerA, basisA, basisB, envDmat, 1.0);

  SPMatrix<SPIN> f_AAA1 = abExchangePot.getMatrix();

  EXPECT_NEAR(f_AAA1(0, 0), -0.55778539409194317, 1e-5);
  EXPECT_NEAR(f_AAA1(1, 0), -0.45434649785726505, 1e-5);
  EXPECT_NEAR(f_AAA1(2, 0), -0.45643495903195952, 1e-5);
  EXPECT_NEAR(f_AAA1(0, 1), -0.45434649785726511, 1e-5);
  EXPECT_NEAR(f_AAA1(0, 2), -0.45643495903195958, 1e-5);
  EXPECT_NEAR(f_AAA1(0, 3), -0.42905604373708817, 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AB matrix for the density of C.
 */
TEST_F(ABExchangePotentialTest, H2_dimer_uFockMatrixABC) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerC->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABExchangePotential<SPIN> abExchangePot(systemControllerA, basisA, basisB, envDmat, 1.0);

  SPMatrix<SPIN> f_AAA1 = abExchangePot.getMatrix();

  EXPECT_NEAR(f_AAA1.alpha(0, 0), -0.55778539409194317, 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(1, 0), -0.45434649785726505, 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(2, 0), -0.45643495903195952, 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 1), -0.45434649785726511, 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 2), -0.45643495903195958, 1e-5);
  EXPECT_NEAR(f_AAA1.alpha(0, 3), -0.42905604373708817, 1e-5);

  EXPECT_NEAR(f_AAA1.beta(0, 0), -0.55778539409194317, 1e-5);
  EXPECT_NEAR(f_AAA1.beta(1, 0), -0.45434649785726505, 1e-5);
  EXPECT_NEAR(f_AAA1.beta(2, 0), -0.45643495903195952, 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 1), -0.45434649785726511, 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 2), -0.45643495903195958, 1e-5);
  EXPECT_NEAR(f_AAA1.beta(0, 3), -0.42905604373708817, 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
