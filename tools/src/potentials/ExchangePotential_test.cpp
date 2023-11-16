/**
 * @file ExchangePotential_test.cpp
 *
 * @author Moritz Bensberg
 * @date Sep 24, 2019
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
#include "potentials/ExchangePotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/ERIPotential.h"
#include "potentials/Potential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ExchangePotentialTest : public ::testing::Test {
 protected:
  ExchangePotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs)) {
  }

  virtual ~ExchangePotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
};

/**
 * @test ExchangePotentialTest
 * @brief Tests the restricted Fock Matrix of Water
 */
TEST_F(ExchangePotentialTest, Water_rFockMatrix) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  ERIPotential<RESTRICTED> referenceCX(systemController, dMat, 1.0,
                                       systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  ERIPotential<RESTRICTED> referenceC(systemController, dMat, 0.0,
                                      systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  FockMatrix<RESTRICTED> f_ref = referenceCX.getMatrix() - referenceC.getMatrix();
  ExchangePotential<RESTRICTED> exchangeOnly(systemController, dMat, 1.0,
                                             systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  FockMatrix<RESTRICTED> f_x = exchangeOnly.getMatrix();
  double maxDiff = (f_ref - f_x).array().abs().maxCoeff();
  EXPECT_NEAR(0.0, maxDiff, 1e-8);
}

/**
 * @test ExchangePotentialTest
 * @brief Tests the unrestricted Fock Matrix of Water
 */
TEST_F(ExchangePotentialTest, Water_uFockMatrix) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  ERIPotential<UNRESTRICTED> referenceCX(systemController, dMat, 1.0,
                                         systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  ERIPotential<UNRESTRICTED> referenceC(systemController, dMat, 0.0,
                                        systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  FockMatrix<UNRESTRICTED> f_ref = referenceCX.getMatrix() - referenceC.getMatrix();
  ExchangePotential<UNRESTRICTED> exchangeOnly(systemController, dMat, 1.0,
                                               systemController->getSettings().basis.integralThreshold, 0.0, 0.0, 0);
  FockMatrix<UNRESTRICTED> f_x = exchangeOnly.getMatrix();
  FockMatrix<UNRESTRICTED> f_diff = f_ref - f_x;
  double maxDiff_alpha = f_diff.alpha.array().abs().maxCoeff();
  double maxDiff_beta = f_diff.beta.array().abs().maxCoeff();
  EXPECT_NEAR(0.0, maxDiff_alpha, 1e-8);
  EXPECT_NEAR(0.0, maxDiff_beta, 1e-8);
}

TEST_F(ExchangePotentialTest, gradients_restricted) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  ERIPotential<Options::SCF_MODES::RESTRICTED> hf2xPot(systemController, dMatController, 2.0, 1.0e-20, 0.0, 0.0, 0);
  ERIPotential<Options::SCF_MODES::RESTRICTED> hf1xPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0);
  ExchangePotential<Options::SCF_MODES::RESTRICTED> xPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0);
  auto ref = (hf2xPot.getGeomGradients() - hf1xPot.getGeomGradients()).eval();
  auto grad = xPot.getGeomGradients();
  ASSERT_NEAR((ref - grad).cwiseAbs().sum(), 0.0, 1e-8);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

TEST_F(ExchangePotentialTest, gradients_unrestricted) {
  // ToDo: Test open-shell system...
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  ERIPotential<Options::SCF_MODES::UNRESTRICTED> hf2xPot(systemController, dMatController, 2.0, 1.0e-20, 0.0, 0.0, 0);
  ERIPotential<Options::SCF_MODES::UNRESTRICTED> hf1xPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0);
  ExchangePotential<Options::SCF_MODES::UNRESTRICTED> xPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0);
  auto ref = (hf2xPot.getGeomGradients() - hf1xPot.getGeomGradients()).eval();
  auto grad = xPot.getGeomGradients();
  ASSERT_NEAR((ref - grad).cwiseAbs().sum(), 0.0, 1e-8);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

} /* namespace Serenity */
