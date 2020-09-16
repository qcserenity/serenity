/**
 * @file ABCoreHamiltonian_test.cpp
 *
 * @date Mai 15, 2018
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
#include "potentials/ABFockMatrixConstruction/ABCoreHamiltonian.h"
#include "potentials/HCorePotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABCoreHamiltonianTest : public ::testing::Test {
 protected:
  ABCoreHamiltonianTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~ABCoreHamiltonianTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
};

/**
 * @test
 * @brief Tests the construction of the AA matrix. This is idenitcal to what the HCorePotential does.
 */
TEST_F(ABCoreHamiltonianTest, H2_rFockMatrix) {
  const auto& basisA = systemControllerA->getBasisController();
  ABCoreHamiltonian<Options::SCF_MODES::RESTRICTED> abHCore(basisA, basisA, systemControllerA->getGeometry());
  SPMatrix<Options::SCF_MODES::RESTRICTED> f_AA_ABCore = abHCore.getMatrix();
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCoreAA(systemControllerA);
  SPMatrix<Options::SCF_MODES::RESTRICTED> f_AACore = hCoreAA.getMatrix();
  auto diff = (f_AA_ABCore - f_AACore).array().abs().maxCoeff();
  // Allow only numerical noise!
  EXPECT_NEAR(diff, 0.0000000000, 1e-15);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test
 * @brief Tests the construction of the AB matrix.
 */
TEST_F(ABCoreHamiltonianTest, H2_rFockMatrixAB) {
  const auto& basisA = systemControllerA->getBasisController();
  const auto& basisB = systemControllerB->getBasisController();
  ABCoreHamiltonian<Options::SCF_MODES::RESTRICTED> abHCore(basisA, basisB, systemControllerA->getGeometry());
  SPMatrix<Options::SCF_MODES::RESTRICTED> f = abHCore.getMatrix();

  EXPECT_NEAR(f(0, 0), -2.121731152335e-01, 1e-5);
  EXPECT_NEAR(f(1, 0), -3.833260245440e-01, 1e-5);
  EXPECT_NEAR(f(2, 0), -3.704729956227e-01, 1e-5);
  EXPECT_NEAR(f(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(f(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(f(5, 0), -4.110640101472e-01, 1e-5);
  EXPECT_NEAR(f(6, 0), -6.592212408032e-01, 1e-5);
  EXPECT_NEAR(f(7, 0), -6.121522196064e-01, 1e-5);
  EXPECT_NEAR(f(8, 0), -4.889753319210e-01, 1e-5);
  EXPECT_NEAR(f(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(f(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(f(11, 0), 2.877198089046e-01, 1e-5);
  EXPECT_NEAR(f(1, 1), -6.212570746688e-01, 1e-5);
  EXPECT_NEAR(f(2, 2), -2.073187362950e-01, 1e-5);
  EXPECT_NEAR(f(3, 3), 0.0, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

/**
 * @test
 * @brief Tests the construction of the AB matrix.
 */
TEST_F(ABCoreHamiltonianTest, H2_uFockMatrixAB) {
  const auto& basisA = systemControllerA->getBasisController();
  const auto& basisB = systemControllerB->getBasisController();
  ABCoreHamiltonian<Options::SCF_MODES::UNRESTRICTED> abHCore(basisA, basisB, systemControllerA->getGeometry());
  SPMatrix<Options::SCF_MODES::UNRESTRICTED> f = abHCore.getMatrix();

  EXPECT_NEAR(f.alpha(0, 0), -2.121731152335e-01, 1e-5);
  EXPECT_NEAR(f.alpha(1, 0), -3.833260245440e-01, 1e-5);
  EXPECT_NEAR(f.alpha(2, 0), -3.704729956227e-01, 1e-5);
  EXPECT_NEAR(f.alpha(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.alpha(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.alpha(5, 0), -4.110640101472e-01, 1e-5);
  EXPECT_NEAR(f.alpha(6, 0), -6.592212408032e-01, 1e-5);
  EXPECT_NEAR(f.alpha(7, 0), -6.121522196064e-01, 1e-5);
  EXPECT_NEAR(f.alpha(8, 0), -4.889753319210e-01, 1e-5);
  EXPECT_NEAR(f.alpha(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.alpha(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.alpha(11, 0), 2.877198089046e-01, 1e-5);
  EXPECT_NEAR(f.alpha(1, 1), -6.212570746688e-01, 1e-5);
  EXPECT_NEAR(f.alpha(2, 2), -2.073187362950e-01, 1e-5);
  EXPECT_NEAR(f.alpha(3, 3), 0.0, 1e-5);

  EXPECT_NEAR(f.beta(0, 0), -2.121731152335e-01, 1e-5);
  EXPECT_NEAR(f.beta(1, 0), -3.833260245440e-01, 1e-5);
  EXPECT_NEAR(f.beta(2, 0), -3.704729956227e-01, 1e-5);
  EXPECT_NEAR(f.beta(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.beta(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.beta(5, 0), -4.110640101472e-01, 1e-5);
  EXPECT_NEAR(f.beta(6, 0), -6.592212408032e-01, 1e-5);
  EXPECT_NEAR(f.beta(7, 0), -6.121522196064e-01, 1e-5);
  EXPECT_NEAR(f.beta(8, 0), -4.889753319210e-01, 1e-5);
  EXPECT_NEAR(f.beta(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.beta(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(f.beta(11, 0), 2.877198089046e-01, 1e-5);
  EXPECT_NEAR(f.beta(1, 1), -6.212570746688e-01, 1e-5);
  EXPECT_NEAR(f.beta(2, 2), -2.073187362950e-01, 1e-5);
  EXPECT_NEAR(f.beta(3, 3), 0.0, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

} /* namespace Serenity */
