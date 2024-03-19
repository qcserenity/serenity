/**
 * @file ABLRExchangePotential_test.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 9, 2020
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
#include "potentials/ABFockMatrixConstruction/ABLRExchangePotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "potentials/LRXPotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABLRExchangePotentialTest : public ::testing::Test {
 protected:
  ABLRExchangePotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)) {
  }

  virtual ~ABLRExchangePotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
};

/**
 * @test
 * @brief Tests the restricted construction of the f_AA matrix for the density of A.
 */
TEST_F(ABLRExchangePotentialTest, H2_dimer_rFockMatrixAAA) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  //  auto basisB = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerA->getElectronicStructure<SPIN>()->getDensityMatrixController());

  ABLRExchangePotential<SPIN> abLRExchangePot(systemControllerA, basisA, basisA, envDmat, 0.7, 0.33);

  LRXPotential<SPIN> lrxPot(systemControllerA, envDmat[0], 0.7, 1e-10, 1e-10, 1e-10, 0, 0.33);

  SPMatrix<SPIN> f_AAA1 = abLRExchangePot.getMatrix();
  SPMatrix<SPIN> f_AAA2 = lrxPot.getMatrix();
  SPMatrix<SPIN> test = f_AAA1.transpose();
  // Test symmetry.
  EXPECT_NEAR((f_AAA1 - test).array().abs().maxCoeff(), 0.0, 1e-12);
  // Test entries.
  double diff = (f_AAA1 - f_AAA2).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-8);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
