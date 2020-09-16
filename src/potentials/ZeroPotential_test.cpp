/**
 * @file ZeroPotential_test.cpp
 *
 * @date Dec 14, 2016
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
#include "potentials/ZeroPotential.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/Potential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ZeroPotentialTest : public ::testing::Test {
 protected:
  ZeroPotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~ZeroPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
};

/**
 * @test ZeroPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(ZeroPotentialTest, H2_dimer_rFockMatrix) {
  auto basisController = systemController->getAtomCenteredBasisController();
  ZeroPotential<Options::SCF_MODES::RESTRICTED> zeroPot(basisController);
  FockMatrix<Options::SCF_MODES::RESTRICTED> F(std::move(zeroPot.getMatrix()));

  // No, this is not a joke. It's called completion. Deal with it.
  EXPECT_NEAR(F(0, 0), 0.0, 1e-9);
  EXPECT_NEAR(F(1, 0), 0.0, 1e-9);
  EXPECT_NEAR(F(2, 0), 0.0, 1e-9);
  EXPECT_NEAR(F(0, 1), 0.0, 1e-9);
  EXPECT_NEAR(F(0, 2), 0.0, 1e-9);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test ZeroPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(ZeroPotentialTest, H2_dimer_uFockMatrix) {
  auto basisController = systemController->getAtomCenteredBasisController();
  ZeroPotential<Options::SCF_MODES::UNRESTRICTED> zeroPot(basisController);
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(zeroPot.getMatrix()));

  // No, this is not a joke. It's called completion. Deal with it.
  EXPECT_NEAR(F.alpha(0, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.alpha(1, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.alpha(2, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.alpha(0, 1), 0.0, 1e-9);
  EXPECT_NEAR(F.alpha(0, 2), 0.0, 1e-9);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(0, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(1, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(2, 0), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(0, 1), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(0, 2), 0.0, 1e-9);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test ZeroPotentialTest
 * @brief Tests the energy.
 */
TEST_F(ZeroPotentialTest, H2_dimer_Energy) {
  auto basisController = systemController->getAtomCenteredBasisController();
  ZeroPotential<Options::SCF_MODES::RESTRICTED> zeroPot(basisController);
  DensityMatrix<Options::SCF_MODES::RESTRICTED> p(basisController);
  p(0, 0) = 1.0;
  p(1, 1) = 1.0;
  double e = zeroPot.getEnergy(p);
  EXPECT_NEAR(e, 0.0, 1e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test ZeroPotentialTest
 * @brief Tests the gradients.
 */
TEST_F(ZeroPotentialTest, H2_dimer_Gradients) {
  auto basisController = systemController->getAtomCenteredBasisController();
  ZeroPotential<Options::SCF_MODES::RESTRICTED> zeroPot(basisController);
  auto grad = zeroPot.getGeomGradients();
  EXPECT_NEAR(grad(0, 0), 0.0, 1e-9);
  EXPECT_NEAR(grad(1, 0), 0.0, 1e-9);
  EXPECT_NEAR(grad(0, 1), 0.0, 1e-9);
  EXPECT_NEAR(grad(1, 1), 0.0, 1e-9);
  EXPECT_NEAR(grad(0, 2), 0.0, 1e-9);
  EXPECT_NEAR(grad(1, 2), 0.0, 1e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test PotentialTest
 * @brief Tests the potential basis comparison.
 */
TEST_F(ZeroPotentialTest, CompareBasis) {
  auto basisController = systemController->getAtomCenteredBasisController();
  ZeroPotential<Options::SCF_MODES::RESTRICTED> zeroPot(basisController);
  auto grad = zeroPot.getGeomGradients();
  EXPECT_TRUE(zeroPot.compareBasis(basisController));
  EXPECT_FALSE(zeroPot.compareBasis(nullptr));
  EXPECT_FALSE(zeroPot.compareBasis(nullptr));
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

} // namespace Serenity
