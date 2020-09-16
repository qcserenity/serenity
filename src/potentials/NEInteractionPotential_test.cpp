/**
 * @file NEInteractionPotential_test.cpp
 *
 * @date Dec 5, 2016
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
#include "potentials/NEInteractionPotential.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/matrices/FockMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/wrappers/Libint.h"
#include "potentials/Potential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class NEInteractionPotentialTest : public ::testing::Test {
 protected:
  NEInteractionPotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
      systemControllerAct(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerEnv(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~NEInteractionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
  std::shared_ptr<SystemController> systemControllerAct;
  std::shared_ptr<SystemController> systemControllerEnv;
};

/**
 * @test NEInteractionPotentialTest
 * @brief Tests the Fock Matrix of a ficticious H2 dimer
 */
TEST_F(NEInteractionPotentialTest, H2_dimer_rFockMatrix) {
  std::shared_ptr<BasisController> basis = systemController->getAtomCenteredBasisController();
  std::vector<std::shared_ptr<const Geometry>> geometries;

  geometries.push_back(systemController->getGeometry());

  NEInteractionPotential<Options::SCF_MODES::RESTRICTED> NEPot(systemController, {systemController}, basis, geometries);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = NEPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), -2.82884, 1e-5);
  EXPECT_NEAR(F(1, 0), -1.56022, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.775952, 1e-5);
  EXPECT_NEAR(F(0, 1), -1.56022, 1e-5);
  EXPECT_NEAR(F(0, 2), -0.775952, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);
}

/**
 * @test NEInteractionPotentialTest
 * @brief Tests the Fock Matrix of a ficticious H2 dimer
 */
TEST_F(NEInteractionPotentialTest, H2_dimer_uFockMatrix) {
  std::shared_ptr<BasisController> basis = systemController->getAtomCenteredBasisController();
  std::vector<std::shared_ptr<const Geometry>> geometries;

  geometries.push_back(systemController->getGeometry());

  NEInteractionPotential<Options::SCF_MODES::UNRESTRICTED> NEPot(systemController, {systemController}, basis, geometries);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = std::move(NEPot.getMatrix());

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), -2.82884, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), -1.56022, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), -0.775952, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), -1.56022, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), -0.775952, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), -2.82884, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), -1.56022, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), -0.775952, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), -1.56022, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), -0.775952, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);
}

/**
 * @test NEInteractionPotentialTest
 * @brief Tests the Fock Matrix of a fictitious H2 dimer
 */
TEST_F(NEInteractionPotentialTest, H2_dimer_rgrad) {
  std::shared_ptr<BasisController> basisAct = systemControllerAct->getAtomCenteredBasisController();
  std::vector<std::shared_ptr<const Geometry>> geometries;

  geometries.push_back(systemControllerEnv->getGeometry());

  NEInteractionPotential<Options::SCF_MODES::RESTRICTED> NEPot(systemControllerAct, {systemControllerEnv}, basisAct, geometries);

  auto NucEl = NEPot.getGeomGradients();

  EXPECT_NEAR(NucEl(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(0, 2), -1.0113744689574913, 1e-5);
  EXPECT_NEAR(NucEl(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(1, 2), -0.97097493046364303, 1e-5);
}

/**
 * @test NEInteractionPotentialTest
 * @brief Tests the Fock Matrix of a ficticious H2 dimer
 */
TEST_F(NEInteractionPotentialTest, H2_dimer_ugrad) {
  std::shared_ptr<BasisController> basisAct = systemControllerAct->getAtomCenteredBasisController();
  std::vector<std::shared_ptr<const Geometry>> geometries;

  geometries.push_back(systemControllerEnv->getGeometry());

  NEInteractionPotential<Options::SCF_MODES::UNRESTRICTED> NEPot(systemControllerAct, {systemControllerEnv}, basisAct,
                                                                 geometries);

  auto NucEl = NEPot.getGeomGradients();

  EXPECT_NEAR(NucEl(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(0, 2), -0.96838583691439761, 1e-5);
  EXPECT_NEAR(NucEl(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(NucEl(1, 2), -1.0157081522761979, 1e-5);
}

} // namespace Serenity
