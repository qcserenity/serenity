/**
 * @file ABFuncPotential_test.cpp
 *
 * @date May 16, 2018
 * @author Moritz Bensberg
 *
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
#include "potentials/ABFockMatrixConstruction/ABFuncPotential.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ExternalDensityOnGridController.h"
#include "data/matrices/DensityMatrixController.h"
#include "potentials/FuncPotential.h"
#include "settings/DFTOptions.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {
class ABFuncPotentialTest : public ::testing::Test {
 protected:
  ABFuncPotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~ABFuncPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
};

/**
 * @test ABFuncPotentialTest
 * @brief Tests the AA Fock Matrix of an H2 (LDA). This is identical to the equivalent test in FuncPotential_test.cpp.
 */
TEST_F(ABFuncPotentialTest, H2_FockMatrixAA_LDA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  auto basisA = systemControllerA->getBasisController();
  auto gridAA = systemControllerA->getGridController();

  ABFuncPotential<Options::SCF_MODES::RESTRICTED> abFuncPot(
      systemControllerA, basisA, basisA, gridAA, {dMat},
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(
      systemControllerA, dMat, gridAA, CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F_AB = abFuncPot.getMatrix();
  SPMatrix<Options::SCF_MODES::RESTRICTED> F_AA = funcPot.getMatrix();

  auto diff = (F_AB - F_AA).array().abs().maxCoeff();
  // allow only white noise!
  EXPECT_NEAR(diff, 0.0, 1e-12);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}
/**
 * @test ABFuncPotentialTest
 * @brief Tests the AA Fock Matrix of an H2 (GGA). This is identical to the equivalent test in FuncPotential_test.cpp.
 */
TEST_F(ABFuncPotentialTest, H2_FockMatrixAA_GGA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  auto basisA = systemControllerA->getBasisController();
  auto gridAA = systemControllerA->getGridController();

  ABFuncPotential<Options::SCF_MODES::RESTRICTED> abFuncPot(
      systemControllerA, basisA, basisA, gridAA, {dMat},
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F = abFuncPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), -0.61072924201172829, 1e-3);
  EXPECT_NEAR(F(1, 0), -0.37691230665476133, 1e-3);
  EXPECT_NEAR(F(2, 0), -0.19074153073625871, 1e-3);
  EXPECT_NEAR(F(0, 1), -0.37691230665476133, 1e-3);
  EXPECT_NEAR(F(0, 2), -0.19074153073625871, 1e-3);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-3);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}
/**
 * @test ABFuncPotentialTest
 * @brief Tests the AB Fock Matrix of an H2 (LDA).
 */
TEST_F(ABFuncPotentialTest, H2_FockMatrixAB_LDA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerA->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  auto gridAA = systemControllerA->getGridController();

  ABFuncPotential<Options::SCF_MODES::RESTRICTED> abFuncPot(
      systemControllerA, basisA, basisB, gridAA, {dMat},
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA));

  SPMatrix<Options::SCF_MODES::RESTRICTED> F = abFuncPot.getMatrix();
  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), -3.636117713288e-02, 1e-5);
  EXPECT_NEAR(F(1, 0), -1.001498923041e-01, 1e-5);
  EXPECT_NEAR(F(2, 0), -1.061975191271e-01, 1e-5);
  EXPECT_NEAR(F(0, 1), -1.116016285240e-01, 1e-5);
  EXPECT_NEAR(F(0, 2), -6.955048960966e-03, 1e-5);
  EXPECT_NEAR(F(0, 3), -6.493820327824e-02, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}
/**
 * @test ABFuncPotentialTest
 * @brief Tests the unrestricted AB Fock Matrix of an H2 (GGA).
 */
TEST_F(ABFuncPotentialTest, H2_uFockMatrixAB_GGA) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemControllerA->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  auto gridAA = systemControllerA->getGridController();

  ABFuncPotential<Options::SCF_MODES::UNRESTRICTED> abFuncPot(
      systemControllerA, basisA, basisB, gridAA, {dMat},
      CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::BP86));

  SPMatrix<Options::SCF_MODES::UNRESTRICTED> F = abFuncPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), -3.801470612146e-02, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), -1.037127215935e-01, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), -1.091671374726e-01, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), -1.165296890920e-01, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), -7.273116668792e-03, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), -6.778065332781e-02, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), -3.801470612146e-02, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), -1.037127215935e-01, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), -1.091671374726e-01, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), -1.165296890920e-01, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), -7.273116668792e-03, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), -6.778065332781e-02, 1e-5);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
}

} /* namespace Serenity */
