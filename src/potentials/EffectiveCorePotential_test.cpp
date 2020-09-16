/**
 * @file EffectiveCorePotential_test.cpp
 *
 *  @date Jun 12, 2018
 *  @author Moritz Bensberg
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
#include "potentials/EffectiveCorePotential.h"
#include "data/matrices/FockMatrix.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {

class EffectiveCorePotentialTest : public ::testing::Test {
 protected:
  EffectiveCorePotentialTest()
    : systemController_HI(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE)) {
  }

  virtual ~EffectiveCorePotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController_HI;
};

/**
 * @test EffectiveCorePotentialTest
 * @brief Tests the ECP-Fock Matrix of HI via the one electron integral controller.
 */
TEST_F(EffectiveCorePotentialTest, HI_rECPMatrix) {
  FockMatrix<Options::SCF_MODES::RESTRICTED> F = systemController_HI->getOneElectronIntegralController()->getECPIntegrals();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 3.6134546941301289, 1e-5);
  EXPECT_NEAR(F(1, 0), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F(2, 0), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(5, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(6, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(0, 1), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F(0, 2), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F(0, 3), 1.1233694306815634, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
}

/**
 * @test EffectiveCorePotentialTest
 * @brief Tests the ECP-Fock Matrix of HI via the EffectiveCorePotental.h
 */
TEST_F(EffectiveCorePotentialTest, HI_rECPMatrixEffectiveCorePotentialClass) {
  EffectiveCorePotential<Options::SCF_MODES::RESTRICTED> ecp(
      systemController_HI, systemController_HI->getGeometry()->getAtoms(), systemController_HI->getBasisController());

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = ecp.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 3.6134546941301289, 1e-5);
  EXPECT_NEAR(F(1, 0), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F(2, 0), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(5, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(6, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(0, 1), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F(0, 2), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F(0, 3), 1.1233694306815634, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
}

/**
 * @test EffectiveCorePotentialTest
 * @brief Tests the unrestricted ECP-Fock Matrix of HI via the EffectiveCorePotental.h
 */
TEST_F(EffectiveCorePotentialTest, HI_uECPMatrixEffectiveCorePotentialClass) {
  EffectiveCorePotential<Options::SCF_MODES::UNRESTRICTED> ecp(
      systemController_HI, systemController_HI->getGeometry()->getAtoms(), systemController_HI->getBasisController());

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = ecp.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 3.6134546941301289, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F.alpha(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(5, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(6, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 1.1233694306815634, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 3.6134546941301289, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F.beta(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(5, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(6, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 6.3833005006243999, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 2.1951043848535425, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 1.1233694306815634, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
}

} /* namespace Serenity */
