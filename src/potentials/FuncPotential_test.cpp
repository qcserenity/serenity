/**
 * @file FuncPotential_test.cpp
 * @author: Kevin Klahr
 *
 * @date 1. December 2016
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/FuncPotential.h"
#include "input/FunctionalClassResolver.h"
#include "potentials/Potential.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "dft/functionals/wrappers/XCFun.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class FuncPotentialTest : public ::testing::Test {
protected:
	FuncPotentialTest(){
	}

	virtual ~FuncPotentialTest() = default;

	  static void TearDownTestCase() {
	    SystemController__TEST_SUPPLY::cleanUp();
	  }


};

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_LDA) {

  auto systemController =
              SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = funcPot.getMatrix();

  // TODO create more data for the test
    EXPECT_NEAR(F(0,0), -0.57441271, 1e-3);
    EXPECT_NEAR(F(1,0), -0.35994699, 1e-3);
    EXPECT_NEAR(F(2,0), -0.18209255, 1e-3);
    EXPECT_NEAR(F(0,1), -0.35994699, 1e-3);
    EXPECT_NEAR(F(0,2), -0.18209255, 1e-3);
    EXPECT_NEAR(F(0,3), 0.0, 1e-3);

    SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_LDA_UNRES) {

  auto systemController =
              SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::LDA));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(funcPot.getMatrix()));

  // TODO create more data for the test
    EXPECT_NEAR(F.alpha(0,0), -0.57441271, 1e-3);
    EXPECT_NEAR(F.alpha(1,0), -0.35994699, 1e-3);
    EXPECT_NEAR(F.alpha(2,0), -0.18209255, 1e-3);
    EXPECT_NEAR(F.alpha(0,1), -0.35994699, 1e-3);
    EXPECT_NEAR(F.alpha(0,2), -0.18209255, 1e-3);
    EXPECT_NEAR(F.alpha(0,3), 0.0, 1e-3);
    EXPECT_NEAR(F.beta(0,0), -0.57441271, 1e-3);
    EXPECT_NEAR(F.beta(1,0), -0.35994699, 1e-3);
    EXPECT_NEAR(F.beta(2,0), -0.18209255, 1e-3);
    EXPECT_NEAR(F.beta(0,1), -0.35994699, 1e-3);
    EXPECT_NEAR(F.beta(0,2), -0.18209255, 1e-3);
    EXPECT_NEAR(F.beta(0,3), 0.0, 1e-3);

    SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_GGA) {

  auto systemController =
              SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = funcPot.getMatrix();

  // TODO create more data for the test
    EXPECT_NEAR(F(0,0), -0.61072924201172829, 1e-3);
    EXPECT_NEAR(F(1,0), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F(2,0), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F(0,1), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F(0,2), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F(0,3), 0.0, 1e-3);

    SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the Fock Matrix of an unrestricted H2 dimer
 */
TEST_F(FuncPotentialTest, H2_FockMatrix_GGA_UNRES) {

  auto systemController =
              SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86));

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(funcPot.getMatrix()));

  // TODO create more data for the test
    EXPECT_NEAR(F.alpha(0,0), -0.61072924201172829, 1e-3);
    EXPECT_NEAR(F.alpha(1,0), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F.alpha(2,0), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F.alpha(0,1), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F.alpha(0,2), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F.alpha(0,3), 0.0, 1e-3);
    EXPECT_NEAR(F.beta(0,0), -0.61072924201172829, 1e-3);
    EXPECT_NEAR(F.beta(1,0), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F.beta(2,0), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F.beta(0,1), -0.37691230665476133, 1e-3);
    EXPECT_NEAR(F.beta(0,2), -0.19074153073625871, 1e-3);
    EXPECT_NEAR(F.beta(0,3), 0.0, 1e-3);

    SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the LDA XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_LDA) {

  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::LDA));
  auto result = funcPot.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0,0),  0.0, 1e-5);
  EXPECT_NEAR(result(0,1),  0.0, 1e-5);
  EXPECT_NEAR(result(0,2), -0.24575394748581764, 1e-4);
  EXPECT_NEAR(result(1,0),  0.0, 1e-5);
  EXPECT_NEAR(result(1,1),  0.0, 1e-5);
  EXPECT_NEAR(result(1,2),  0.24575394748581764, 1e-4);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);
}

/**
 * @test FuncPotentialTest
 * @brief Tests the unrestricted LDA XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_LDA_UNRES) {

  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::LDA));
  auto result = funcPot.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0,0),  0.0, 1e-5);
  EXPECT_NEAR(result(0,1),  0.0, 1e-5);
  EXPECT_NEAR(result(0,2), -0.24575394748581764, 1e-4);
  EXPECT_NEAR(result(1,0),  0.0, 1e-5);
  EXPECT_NEAR(result(1,1),  0.0, 1e-5);
  EXPECT_NEAR(result(1,2),  0.24575394748581764, 1e-4);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the BP86 XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_GGA) {

  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  FuncPotential<Options::SCF_MODES::RESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86));
  auto result = funcPot.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0,0), 0.0, 1e-5);
  EXPECT_NEAR(result(0,1), 0.0, 1e-5);
  EXPECT_NEAR(result(0,2), -0.2542114, 1e-4);
  EXPECT_NEAR(result(1,0), 0.0, 1e-5);
  EXPECT_NEAR(result(1,1), 0.0, 1e-5);
  EXPECT_NEAR(result(1,2), 0.2542114, 1e-4);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the unrestricted BP86 XC part of the gradient.
 */
TEST_F(FuncPotentialTest, H2_Gradient_GGA_UNRES) {

  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  FuncPotential<Options::SCF_MODES::UNRESTRICTED> funcPot(systemController, dMat, grid, FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86));
  auto result = funcPot.getGeomGradients();

  // Reference for this test are XC Gradients as of 22.03.16
  EXPECT_NEAR(result(0,0), 0.0, 1e-5);
  EXPECT_NEAR(result(0,1), 0.0, 1e-5);
  EXPECT_NEAR(result(0,2), -0.2542114, 1e-4);
  EXPECT_NEAR(result(1,0), 0.0, 1e-5);
  EXPECT_NEAR(result(1,1), 0.0, 1e-5);
  EXPECT_NEAR(result(1,2), 0.2542114, 1e-4);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);
}
/**
 * @test FuncPotentialTest
 * @brief Tests the BP86 matrix-potential calculated 2 ways.
 */
TEST_F(FuncPotentialTest, H2_Potential_GGA) {

  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();

  std::shared_ptr<GridController> grid = systemController->getGridController();

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), dMat->getDensityMatrix().getBasisController(), grid);
  basisFunctionOnGridController->setHighestDerivative(2);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED> >(
      basisFunctionOnGridController, 0.0);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED> >(
      densOnGridCalculator, dMat, 2);
  auto gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED> >(
      basisFunctionOnGridController, 0.0);

  XCFun<Options::SCF_MODES::RESTRICTED> xcFun(128);
  auto funcData1 = xcFun.calcData(
      FUNCTIONAL_DATA_TYPE::GRADIENTS,
      FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86),
      densOnGridController);

  auto funcData2 = xcFun.calcData(
      FUNCTIONAL_DATA_TYPE::POTENTIAL,
      FunctionalClassResolver::resolveFunctional(Options::FUNCTIONALS::BP86),
      densOnGridController);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F1(dMat->getDensityMatrix().getBasisController());
  FockMatrix<Options::SCF_MODES::RESTRICTED> F2(dMat->getDensityMatrix().getBasisController());

  gridToMatrix->addScalarOperatorToMatrix(F2,*funcData2.potential);
  gridToMatrix->addScalarOperatorToMatrix(F1,*funcData1.dFdRho,*funcData1.dFdGradRho);

  EXPECT_NEAR(F1(0,0), F2(0,0), 3e-6);
  EXPECT_NEAR(F1(1,0), F2(1,0), 1e-6);
  EXPECT_NEAR(F1(1,1), F2(1,1), 1e-6);
  EXPECT_NEAR(F1(2,0), F2(2,0), 3e-6);
  EXPECT_NEAR(F1(2,1), F2(2,1), 1e-6);
  EXPECT_NEAR(F1(2,2), F2(2,2), 3e-6);
  EXPECT_NEAR(F1(3,0), F2(3,0), 1e-6);
  EXPECT_NEAR(F1(3,1), F2(3,1), 1e-6);
  EXPECT_NEAR(F1(3,2), F2(3,2), 1e-6);
  EXPECT_NEAR(F1(3,3), F2(3,3), 1e-6);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86);
}
} /*Namespace*/
