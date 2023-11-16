/**
 * @file LibXC_test.cpp
 * @author: Jan Unsleber
 *
 * @date Sep 20, 2017
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
#ifdef SERENITY_USE_LIBXC
/* Include Serenity Internal Headers */
#include "dft/functionals/wrappers/LibXC.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class LibXCTest : public ::testing::Test {
 protected:
  LibXCTest() {
  }
  virtual ~LibXCTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

#if defined SERENITY_USE_LIBXC && defined SERENITY_USE_XCFUN
/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_LDA_Restricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  LibXC<RESTRICTED> libxc(128);
  XCFun<RESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<RESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<RESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA);

  // Calculate data
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 1e-9);
  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-9);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    ASSERT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*funResults.dFdRho)[i], (*libResults.dFdRho)[i], 1.0e-5);
    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan((*libResults.d2FdRho2)[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_LDA_Unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP);

  LibXC<UNRESTRICTED> libxc(128);
  XCFun<UNRESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<UNRESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA);

  // Calculate data
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 1e-9);
  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-9);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    EXPECT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-5);
    EXPECT_NEAR(funResults.dFdRho->alpha[i], libResults.dFdRho->alpha[i], 1.0e-4);
    EXPECT_NEAR(funResults.dFdRho->beta[i], libResults.dFdRho->beta[i], 1.0e-4);
    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->bb[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Restricted_Invariant) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  LibXC<RESTRICTED> libxc(128);
  XCFun<RESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<RESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<RESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 1e-4);

  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-6);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*funResults.dFdRho)[i], (*libResults.dFdRho)[i], 1.0e-5);
    // Derivatives of sigma are fairly different between libxc and xcfun
    ASSERT_FALSE(std::isnan((*libResults.dFdSigma)[i]));
    // ASSERT_NEAR((*funResults.dFdSigma)[i], (*libResults.dFdSigma)[i], 1.0e-5);
    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan((*libResults.d2FdRho2)[i]));
    ASSERT_FALSE(std::isnan((*libResults.d2FdRhodSigma)[i]));
    ASSERT_FALSE(std::isnan((*libResults.d2FdSigma2)[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Unrestricted_Invariant) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP);

  LibXC<UNRESTRICTED> libxc(128);
  XCFun<UNRESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<UNRESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 0.2);
  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-5);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-3);
    ASSERT_NEAR(funResults.dFdRho->alpha[i], libResults.dFdRho->alpha[i], 1.0e-4);
    ASSERT_NEAR(funResults.dFdRho->beta[i], libResults.dFdRho->beta[i], 1.0e-4);
    // Derivatives of sigma are fairly different between libxc and xcfun
    ASSERT_FALSE(std::isnan(libResults.dFdSigma->aa[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdSigma->ab[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdSigma->bb[i]));
    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->aaa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->aab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->abb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->baa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->bab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodSigma->bbb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->aaaa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->aaab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->aabb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->abab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->abbb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdSigma2->bbbb[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Restricted_Cartesian) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  LibXC<RESTRICTED> libxc(128);
  XCFun<RESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<RESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<RESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 1e-4);
  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-6);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*funResults.dFdRho)[i], (*libResults.dFdRho)[i], 1.0e-5);
    // Derivatives of sigma are fairly different between libxc and xcfun
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->x[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->y[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->z[i]));

    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan((*libResults.d2FdRho2)[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->x[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->y[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->z[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xx[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xy[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xz[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yy[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yz[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->zz[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Unrestricted_Cartesian) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP);

  LibXC<UNRESTRICTED> libxc(128);
  XCFun<UNRESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<UNRESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  // Compare
  ASSERT_NEAR(funResults.epuv->sum(), libResults.epuv->sum(), 0.2);
  ASSERT_NEAR(funResults.energy, libResults.energy, 1e-5);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*funResults.epuv)[i], (*libResults.epuv)[i], 1.0e-3);
    ASSERT_NEAR(funResults.dFdRho->alpha[i], libResults.dFdRho->alpha[i], 1.0e-4);
    ASSERT_NEAR(funResults.dFdRho->beta[i], libResults.dFdRho->beta[i], 1.0e-4);
    // Derivatives of sigma are fairly different between libxc and xcfun
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->x.alpha[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->x.beta[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->y.alpha[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->y.beta[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->z.alpha[i]));
    ASSERT_FALSE(std::isnan(libResults.dFdGradRho->z.beta[i]));

    // 2nd derivatives can be way off for small densities just make sure they are non nan.
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRho2->bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->x.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->y.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->z.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->x.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->y.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->z.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->x.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->y.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->z.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->x.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->y.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdRhodGradRho->z.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xx.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xy.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xz.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xx.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xy.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xz.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yy.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yz.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xy.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yy.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yz.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->zz.aa[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xz.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yz.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->zz.ab[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xx.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xy.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xz.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yy.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yz.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->zz.bb[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->xx.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->yy.ba[i]));
    ASSERT_FALSE(std::isnan(libResults.d2FdGradRho2->zz.ba[i]));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Restricted_Cartesian_Cross) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP);

  LibXC<RESTRICTED> libxc(128);
  XCFun<RESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<RESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<RESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  auto grid = densOnGridController->getGridController();
  libResults.dFdGradRho = makeGradientPtr<GridPotential<RESTRICTED>>(grid);
  libResults.d2FdRhodGradRho = makeGradientPtr<DoublySpinPolarizedData<RESTRICTED, GridData<RESTRICTED>>>(grid);
  libResults.d2FdGradRho2 = makeHessianPtr<DoublySpinPolarizedData<RESTRICTED, GridData<RESTRICTED>>>(grid);
  auto gradient = std::make_shared<Gradient<DensityOnGrid<RESTRICTED>>>(densOnGridController->getDensityGradientOnGrid());
  libxc.complete(libResults, *gradient, 0, grid->getNGridPoints());

  // Compare
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // 1st deriv
    ASSERT_NEAR(funResults.dFdGradRho->x[i], libResults.dFdGradRho->x[i],
                std::max(fabs(funResults.dFdGradRho->x[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->y[i], libResults.dFdGradRho->y[i],
                std::max(fabs(funResults.dFdGradRho->y[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->z[i], libResults.dFdGradRho->z[i],
                std::max(fabs(funResults.dFdGradRho->z[i] * 1.0e-6), 1.0e-5));
    // 2nd deriv
    ASSERT_NEAR((*funResults.d2FdRho2)[i], (*libResults.d2FdRho2)[i],
                std::max(fabs((*funResults.d2FdRho2)[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->x[i], libResults.d2FdRhodGradRho->x[i],
                std::max(fabs(funResults.d2FdRhodGradRho->x[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->y[i], libResults.d2FdRhodGradRho->y[i],
                std::max(fabs(funResults.d2FdRhodGradRho->y[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->z[i], libResults.d2FdRhodGradRho->z[i],
                std::max(fabs(funResults.d2FdRhodGradRho->z[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xx[i], libResults.d2FdGradRho2->xx[i],
                std::max(fabs(funResults.d2FdGradRho2->xx[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xy[i], libResults.d2FdGradRho2->xy[i],
                std::max(fabs(funResults.d2FdGradRho2->xy[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xz[i], libResults.d2FdGradRho2->xz[i],
                std::max(fabs(funResults.d2FdGradRho2->xz[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yy[i], libResults.d2FdGradRho2->yy[i],
                std::max(fabs(funResults.d2FdGradRho2->yy[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yz[i], libResults.d2FdGradRho2->yz[i],
                std::max(fabs(funResults.d2FdGradRho2->yz[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->zz[i], libResults.d2FdGradRho2->zz[i],
                std::max(fabs(funResults.d2FdGradRho2->zz[i] * 1.0e-6), 1.0e-5));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Unrestricted_Cartesian_Cross) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP);

  LibXC<UNRESTRICTED> libxc(128);
  XCFun<UNRESTRICTED> xcfun(128);

  // Prep input
  auto bfOnGrid = BasisFunctionOnGridControllerFactory::produce(act->getSettings(), act->getBasisController(),
                                                                act->getGridController());
  auto densOnGridCalc =
      std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(bfOnGrid, act->getSettings().grid.blockAveThreshold);
  auto es = act->getElectronicStructure<UNRESTRICTED>();
  auto dMatController = es->getDensityMatrixController();
  auto densOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densOnGridCalc, dMatController);
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto libResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 2);
  auto funResults = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 2);

  auto grid = densOnGridController->getGridController();
  libResults.dFdGradRho = makeGradientPtr<GridPotential<UNRESTRICTED>>(grid);
  libResults.d2FdRhodGradRho = makeGradientPtr<DoublySpinPolarizedData<UNRESTRICTED, GridData<RESTRICTED>>>(grid);
  libResults.d2FdGradRho2 = makeHessianPtr<DoublySpinPolarizedData<UNRESTRICTED, GridData<RESTRICTED>>>(grid);
  auto gradient = std::make_shared<Gradient<DensityOnGrid<UNRESTRICTED>>>(densOnGridController->getDensityGradientOnGrid());
  libxc.complete(libResults, *gradient, 0, grid->getNGridPoints());

  // Compare
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // 1st deriv
    ASSERT_NEAR(funResults.dFdGradRho->x.alpha[i], libResults.dFdGradRho->x.alpha[i],
                std::max(fabs(funResults.dFdGradRho->x.alpha[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->y.alpha[i], libResults.dFdGradRho->y.alpha[i],
                std::max(fabs(funResults.dFdGradRho->y.alpha[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->z.alpha[i], libResults.dFdGradRho->z.alpha[i],
                std::max(fabs(funResults.dFdGradRho->z.alpha[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->x.beta[i], libResults.dFdGradRho->x.beta[i],
                std::max(fabs(funResults.dFdGradRho->x.beta[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->y.beta[i], libResults.dFdGradRho->y.beta[i],
                std::max(fabs(funResults.dFdGradRho->y.beta[i] * 1.0e-6), 1.0e-5));
    ASSERT_NEAR(funResults.dFdGradRho->z.beta[i], libResults.dFdGradRho->z.beta[i],
                std::max(fabs(funResults.dFdGradRho->z.beta[i] * 1.0e-6), 1.0e-5));
    // 2nd deriv (ab cross terms seem a bit more noisy than aa and bb terms)
    ASSERT_NEAR(funResults.d2FdRho2->aa[i], libResults.d2FdRho2->aa[i],
                std::max(fabs(funResults.d2FdRho2->aa[i]) * 1.0e-6, 5.0e-5));
    ASSERT_NEAR(funResults.d2FdRho2->ab[i], libResults.d2FdRho2->ab[i],
                std::max(fabs(funResults.d2FdRho2->ab[i]) * 1.0e-6, 5.0e-5));
    ASSERT_NEAR(funResults.d2FdRho2->bb[i], libResults.d2FdRho2->bb[i],
                std::max(fabs(funResults.d2FdRho2->bb[i]) * 1.0e-6, 5.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->x.aa[i], libResults.d2FdRhodGradRho->x.aa[i],
                std::max(fabs(funResults.d2FdRhodGradRho->x.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->y.aa[i], libResults.d2FdRhodGradRho->y.aa[i],
                std::max(fabs(funResults.d2FdRhodGradRho->y.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->z.aa[i], libResults.d2FdRhodGradRho->z.aa[i],
                std::max(fabs(funResults.d2FdRhodGradRho->z.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->x.ab[i], libResults.d2FdRhodGradRho->x.ab[i],
                std::max(fabs(funResults.d2FdRhodGradRho->x.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->y.ab[i], libResults.d2FdRhodGradRho->y.ab[i],
                std::max(fabs(funResults.d2FdRhodGradRho->y.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->z.ab[i], libResults.d2FdRhodGradRho->z.ab[i],
                std::max(fabs(funResults.d2FdRhodGradRho->z.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->x.ba[i], libResults.d2FdRhodGradRho->x.ba[i],
                std::max(fabs(funResults.d2FdRhodGradRho->x.ba[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->y.ba[i], libResults.d2FdRhodGradRho->y.ba[i],
                std::max(fabs(funResults.d2FdRhodGradRho->y.ba[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->z.ba[i], libResults.d2FdRhodGradRho->z.ba[i],
                std::max(fabs(funResults.d2FdRhodGradRho->z.ba[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->x.bb[i], libResults.d2FdRhodGradRho->x.bb[i],
                std::max(fabs(funResults.d2FdRhodGradRho->x.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->y.bb[i], libResults.d2FdRhodGradRho->y.bb[i],
                std::max(fabs(funResults.d2FdRhodGradRho->y.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdRhodGradRho->z.bb[i], libResults.d2FdRhodGradRho->z.bb[i],
                std::max(fabs(funResults.d2FdRhodGradRho->z.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xx.aa[i], libResults.d2FdGradRho2->xx.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->xx.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yy.aa[i], libResults.d2FdGradRho2->yy.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->yy.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->zz.aa[i], libResults.d2FdGradRho2->zz.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->zz.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xx.bb[i], libResults.d2FdGradRho2->xx.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->xx.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yy.bb[i], libResults.d2FdGradRho2->yy.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->yy.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->zz.bb[i], libResults.d2FdGradRho2->zz.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->zz.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xx.ab[i], libResults.d2FdGradRho2->xx.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->xx.ab[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->yy.ab[i], libResults.d2FdGradRho2->yy.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->yy.ab[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->zz.ab[i], libResults.d2FdGradRho2->zz.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->zz.ab[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->xx.ba[i], libResults.d2FdGradRho2->xx.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->xx.ba[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->yy.ba[i], libResults.d2FdGradRho2->yy.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->yy.ba[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->zz.ba[i], libResults.d2FdGradRho2->zz.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->zz.ba[i]) * 1.0e-6, 1.0e-3));
    ASSERT_NEAR(funResults.d2FdGradRho2->xz.aa[i], libResults.d2FdGradRho2->xz.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->xz.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xy.aa[i], libResults.d2FdGradRho2->xy.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->xy.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yz.aa[i], libResults.d2FdGradRho2->yz.aa[i],
                std::max(fabs(funResults.d2FdGradRho2->yz.aa[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xy.bb[i], libResults.d2FdGradRho2->xy.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->xy.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->yz.bb[i], libResults.d2FdGradRho2->yz.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->yz.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xz.bb[i], libResults.d2FdGradRho2->xz.bb[i],
                std::max(fabs(funResults.d2FdGradRho2->xz.bb[i]) * 1.0e-6, 1.0e-5));
    ASSERT_NEAR(funResults.d2FdGradRho2->xy.ab[i], libResults.d2FdGradRho2->xy.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->xy.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdGradRho2->xz.ab[i], libResults.d2FdGradRho2->xz.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->xz.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdGradRho2->yz.ab[i], libResults.d2FdGradRho2->yz.ab[i],
                std::max(fabs(funResults.d2FdGradRho2->yz.ab[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdGradRho2->xy.ba[i], libResults.d2FdGradRho2->xy.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->xy.ba[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdGradRho2->xz.ba[i], libResults.d2FdGradRho2->xz.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->xz.ba[i]) * 1.0e-6, 1.0e-4));
    ASSERT_NEAR(funResults.d2FdGradRho2->yz.ba[i], libResults.d2FdGradRho2->yz.ba[i],
                std::max(fabs(funResults.d2FdGradRho2->yz.ba[i]) * 1.0e-6, 1.0e-4));
  }
  SystemController__TEST_SUPPLY::cleanUp();
}
#endif /* defined SERENITY_USE_LIBXC && defined SERENITY_USE_XCFUN */
/**
 * @test LibXCTest
 * @brief Tests the functional aliases of the implemented functionals.
 */
// TEST_F(LibXCTest, FunctionalAliases) {
//   LibXC<RESTRICTED> LibXC(1);
//   for (unsigned int i=0;i<static_cast<unsigned int>(BasicFunctionals::BASIC_FUNCTIONALS::SAOP);i++){
//     xc_functional func = xc_new_functional();
//     int iErr = xc_set(func,LibXC.getAlias(static_cast<BasicFunctionals::BASIC_FUNCTIONALS>(i)),1.0);
//     EXPECT_TRUE(iErr==0);
//   }
// }

} /* namespace Serenity */
#endif /* SERENITY_USE_LIBXC */
