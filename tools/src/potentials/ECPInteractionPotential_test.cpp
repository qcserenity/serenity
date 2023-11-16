/**
 * @file ECPInteractionPotential_test.cpp
 *
 * @date Jun 12, 2018
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
#include "potentials/ECPInteractionPotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "geometry/Geometry.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {

/**
 * @class ECPInteractionPotentialTest ECPInteractionPotential_test.cpp
 * @brief Tests ECPInteractionPotential.h/.cpp. Note: This is tested indirectly in
 *        the TDEmbeddingTask_test, too. See TDEmbeddingTaskTest.embeddingMode_Levelshift_ECPs
 */
class ECPInteractionPotentialTest : public ::testing::Test {
 protected:
  ECPInteractionPotentialTest()
    : systemController_HI(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE)),
      systemController_I2(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE)) {
  }

  virtual ~ECPInteractionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController_HI;
  std::shared_ptr<SystemController> systemController_I2;
};

/**
 * @test ECPInteractionPotentialTest
 * @brief Tests the ECP-interaction-Fock Matrix of HI and I2.
 */
TEST_F(ECPInteractionPotentialTest, HI_I2_rECPMatrix) {
  const auto scfmode = Options::SCF_MODES::RESTRICTED;
  ECPInteractionPotential<scfmode> ecpInt(
      systemController_HI, systemController_HI->getGeometry()->getAtoms(), systemController_I2->getGeometry()->getAtoms(),
      {systemController_I2->getElectronicStructure<scfmode>()->getDensityMatrixController()},
      systemController_HI->getBasisController());
  FockMatrix<scfmode> F = ecpInt.getMatrix();
  // The results are obtained with a working implementation.
  // The implementation is verified by comparison to
  // the supersystem results in a top-down calculation.
  // ECPs work with projectors onto spherical harmonics. Thus, a lot of entries are
  // zero due to symmetry.
  EXPECT_NEAR(F(0, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(1, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(2, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(4, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(5, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(6, 0), 0.0, 1e-8);
  EXPECT_NEAR(F(3, 3), 1.197306566606e-07, 1e-8);
  EXPECT_NEAR(F(13, 3), 6.093976035420e-07, 1e-8);
  EXPECT_NEAR(F(14, 3), -3.330690707420e-07, 1e-8);
  EXPECT_NEAR(F(13, 13), 3.113843725064e-06, 1e-8);
  EXPECT_NEAR(F(14, 13), -1.682858260680e-06, 1e-8);

  EXPECT_NEAR(ecpInt.getEnergy(systemController_HI->getElectronicStructure<scfmode>()->getDensityMatrix()),
              7.4618060206549205e-06, 1e-8);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_HI->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_HI->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_I2->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
}

/**
 * @test ECPInteractionPotentialTest
 * @brief Tests the unrestricted ECP-interaction-Fock Matrix of HI and I2.
 */
TEST_F(ECPInteractionPotentialTest, HI_I2_uECPMatrix) {
  const auto scfmode = Options::SCF_MODES::UNRESTRICTED;
  ECPInteractionPotential<scfmode> ecpInt(
      systemController_HI, systemController_HI->getGeometry()->getAtoms(), systemController_I2->getGeometry()->getAtoms(),
      {systemController_I2->getElectronicStructure<scfmode>()->getDensityMatrixController()},
      systemController_HI->getBasisController());

  FockMatrix<scfmode> F = ecpInt.getMatrix();

  // The results are obtained with a working implementation.
  // The implementation is verified by comparison to
  // the supersystem results in a top-down calculation.
  // ECPs work with projectors onto spherical harmonics. Thus, a lot of entries are
  // zero due to symmetry.
  EXPECT_NEAR(F.alpha(0, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(1, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(2, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(4, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(5, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(6, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.alpha(3, 3), 1.197306566606e-07, 1e-8);
  EXPECT_NEAR(F.alpha(13, 3), 6.093976035420e-07, 1e-8);
  EXPECT_NEAR(F.alpha(14, 3), -3.330690707420e-07, 1e-8);
  EXPECT_NEAR(F.alpha(13, 13), 3.113843725064e-06, 1e-8);
  EXPECT_NEAR(F.alpha(14, 13), -1.682858260680e-06, 1e-8);

  EXPECT_NEAR(F.beta(0, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(1, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(2, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(4, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(5, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(6, 0), 0.0, 1e-8);
  EXPECT_NEAR(F.beta(3, 3), 1.197306566606e-07, 1e-8);
  EXPECT_NEAR(F.beta(13, 3), 6.093976035420e-07, 1e-8);
  EXPECT_NEAR(F.beta(14, 3), -3.330690707420e-07, 1e-8);
  EXPECT_NEAR(F.beta(13, 13), 3.113843725064e-06, 1e-8);
  EXPECT_NEAR(F.beta(14, 13), -1.682858260680e-06, 1e-8);

  EXPECT_NEAR(ecpInt.getEnergy(systemController_HI->getElectronicStructure<scfmode>()->getDensityMatrix()),
              7.4618060206549205e-06, 1e-8);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_HI->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_HI->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController_I2->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::HI_Def2_SVP_PBE);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
}

} /* namespace Serenity */
