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
  /**
   * @brief This function handles the fact that especially higher derivatives of the density can take on huge values,
   * especially at gridpoints where the value of the density is small. If the absolute value of the difference between
   * the two values is smaller than the threshold, the assertion has passed. Otherwise, the relative error is asserted
   * against the same threshold. The second case is expected to be relevant for large values. Additionally, no assertion
   * is performed if either r1 or r2 is smaller than screen.
   * @param r1 A density (e.g. alpha)
   * @param r2 Another density (or the same as r1 if the quantity only depends on one spin density)
   * @param v1 First value, to be compared against v2.
   * @param v2 Second value, to be compared against v1.
   * @param screen Screening threshold for the density. For very small values of the density, its derivatives are
   * expected to be numerically noisy and are therefore discarded.
   * @param thresh Threshold either for the absolute deviation of v1 and v2 or for their relative deviation.
   */
  void sAssertNear(const double r1, const double r2, const double& v1, const double& v2, const double screen,
                   const double thresh) {
    if ((r1 > screen) && (r2 > screen)) {
      if (abs(v1 - v2) > thresh) {
        ASSERT_NEAR(v1 / v2, 1.0, thresh);
      }
    }
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
  const DensityOnGrid<RESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA);

  // Calculate data
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 1e-9);
  ASSERT_NEAR(fun.energy, lib.energy, 1e-9);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    ASSERT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*fun.dFdRho)[i], (*lib.dFdRho)[i], 1.0e-5);
    sAssertNear(d[i], d[i], (*lib.d2FdRho2)[i], (*fun.d2FdRho2)[i], 1e-10, 5e-5);
    sAssertNear(d[i], d[i], (*lib.d3FdRho3)[i], (*fun.d3FdRho3)[i], 1e-10, 5e-5);
  }
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
  const DensityOnGrid<UNRESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::LDA);

  // Calculate data
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 1e-9);
  ASSERT_NEAR(fun.energy, lib.energy, 1e-9);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    EXPECT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-5);
    EXPECT_NEAR(fun.dFdRho->alpha[i], lib.dFdRho->alpha[i], 1.0e-4);
    EXPECT_NEAR(fun.dFdRho->beta[i], lib.dFdRho->beta[i], 1.0e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRho2->aa[i], fun.d2FdRho2->aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRho2->ab[i], fun.d2FdRho2->ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRho2->bb[i], fun.d2FdRho2->bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho3->aaa[i], fun.d3FdRho3->aaa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->aab[i], fun.d3FdRho3->aab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->abb[i], fun.d3FdRho3->abb[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho3->bbb[i], fun.d3FdRho3->bbb[i], 1e-10, 5e-4);
  }
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Restricted_Invariant) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "sto-6g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 5;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

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
  const DensityOnGrid<RESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 3);
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 5e-4);

  ASSERT_NEAR(fun.energy, lib.energy, 1e-6);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*fun.dFdRho)[i], (*lib.dFdRho)[i], 5.0e-5);
    sAssertNear(d[i], d[i], (*lib.dFdSigma)[i], (*fun.dFdSigma)[i], 1e-12, 1e-5);
    sAssertNear(d[i], d[i], (*lib.d2FdRho2)[i], (*fun.d2FdRho2)[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], (*lib.d2FdRhodSigma)[i], (*fun.d2FdRhodSigma)[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], (*lib.d2FdSigma2)[i], (*fun.d2FdSigma2)[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], (*lib.d3FdRho3)[i], (*fun.d3FdRho3)[i], 1e-8, 5e-3);
    sAssertNear(d[i], d[i], (*lib.d3FdRho2dSigma)[i], (*fun.d3FdRho2dSigma)[i], 1e-8, 5e-3);
    sAssertNear(d[i], d[i], (*lib.d3FdRhodSigma2)[i], (*fun.d3FdRhodSigma2)[i], 1e-8, 5e-3);
    sAssertNear(d[i], d[i], (*lib.d3FdSigma3)[i], (*fun.d3FdSigma3)[i], 1e-8, 5e-3);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
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
  const DensityOnGrid<UNRESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 3);
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 0.2);
  ASSERT_NEAR(fun.energy, lib.energy, 1e-5);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    ASSERT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-3);
    ASSERT_NEAR(fun.dFdRho->alpha[i], lib.dFdRho->alpha[i], 1.0e-4);
    ASSERT_NEAR(fun.dFdRho->beta[i], lib.dFdRho->beta[i], 1.0e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.dFdSigma->aa[i], fun.dFdSigma->aa[i], 1e-12, 1e-5);
    sAssertNear(d.alpha[i], d.beta[i], lib.dFdSigma->ab[i], fun.dFdSigma->ab[i], 1e-12, 1e-5);
    sAssertNear(d.beta[i], d.beta[i], lib.dFdSigma->bb[i], fun.dFdSigma->bb[i], 1e-12, 1e-5);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRho2->aa[i], fun.d2FdRho2->aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRho2->ab[i], fun.d2FdRho2->ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRho2->bb[i], fun.d2FdRho2->bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRhodSigma->aaa[i], fun.d2FdRhodSigma->aaa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodSigma->aab[i], fun.d2FdRhodSigma->aab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodSigma->abb[i], fun.d2FdRhodSigma->abb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodSigma->baa[i], fun.d2FdRhodSigma->baa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodSigma->bab[i], fun.d2FdRhodSigma->bab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRhodSigma->bbb[i], fun.d2FdRhodSigma->bbb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdSigma2->aaaa[i], fun.d2FdSigma2->aaaa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdSigma2->aaab[i], fun.d2FdSigma2->aaab[i], 1e-10, 5e-4);

    ASSERT_NEAR(lib.d2FdSigma2->aaab[i], lib.d2FdSigma2->abbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdRhodSigma2->aaaab[i], lib.d3FdRhodSigma2->aabbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdRhodSigma2->baaab[i], lib.d3FdRhodSigma2->babbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdSigma3->aaaaab[i], lib.d3FdSigma3->aaabbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdSigma3->aaaaab[i], lib.d3FdSigma3->abbbbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdSigma3->aaaabb[i], lib.d3FdSigma3->aabbbb[i], 1e-12);
    ASSERT_NEAR(lib.d3FdSigma3->aaabab[i], lib.d3FdSigma3->ababbb[i], 1e-12);

    ASSERT_NEAR(fun.d2FdSigma2->aaab[i], fun.d2FdSigma2->abbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdRhodSigma2->aaaab[i], fun.d3FdRhodSigma2->aabbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdRhodSigma2->baaab[i], fun.d3FdRhodSigma2->babbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdSigma3->aaaaab[i], fun.d3FdSigma3->aaabbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdSigma3->aaaaab[i], fun.d3FdSigma3->abbbbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdSigma3->aaaabb[i], fun.d3FdSigma3->aabbbb[i], 1e-12);
    ASSERT_NEAR(fun.d3FdSigma3->aaabab[i], fun.d3FdSigma3->ababbb[i], 1e-12);

    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdSigma2->aabb[i], fun.d2FdSigma2->aabb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdSigma2->abab[i], fun.d2FdSigma2->abab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdSigma2->abbb[i], fun.d2FdSigma2->abbb[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdSigma2->bbbb[i], fun.d2FdSigma2->bbbb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho3->aaa[i], fun.d3FdRho3->aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->aab[i], fun.d3FdRho3->aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->abb[i], fun.d3FdRho3->abb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho3->bbb[i], fun.d3FdRho3->bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho2dSigma->aaaa[i], fun.d3FdRho2dSigma->aaaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->aaab[i], fun.d3FdRho2dSigma->aaab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->aabb[i], fun.d3FdRho2dSigma->aabb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->abaa[i], fun.d3FdRho2dSigma->abaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->abab[i], fun.d3FdRho2dSigma->abab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->abbb[i], fun.d3FdRho2dSigma->abbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->bbaa[i], fun.d3FdRho2dSigma->bbaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dSigma->bbab[i], fun.d3FdRho2dSigma->bbab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho2dSigma->bbbb[i], fun.d3FdRho2dSigma->bbbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodSigma2->aaaaa[i], fun.d3FdRhodSigma2->aaaaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->aaaab[i], fun.d3FdRhodSigma2->aaaab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->aaabb[i], fun.d3FdRhodSigma2->aaabb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->aabab[i], fun.d3FdRhodSigma2->aabab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->aabbb[i], fun.d3FdRhodSigma2->aabbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->abbbb[i], fun.d3FdRhodSigma2->abbbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodSigma2->baaaa[i], fun.d3FdRhodSigma2->baaaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->baaab[i], fun.d3FdRhodSigma2->baaab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->baabb[i], fun.d3FdRhodSigma2->baabb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->babab[i], fun.d3FdRhodSigma2->babab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodSigma2->babbb[i], fun.d3FdRhodSigma2->babbb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodSigma2->bbbbb[i], fun.d3FdRhodSigma2->bbbbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdSigma3->aaaaaa[i], fun.d3FdSigma3->aaaaaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->aaaaab[i], fun.d3FdSigma3->aaaaab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->aaaabb[i], fun.d3FdSigma3->aaaabb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->aaabab[i], fun.d3FdSigma3->aaabab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->aaabbb[i], fun.d3FdSigma3->aaabbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->aabbbb[i], fun.d3FdSigma3->aabbbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->ababab[i], fun.d3FdSigma3->ababab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->ababbb[i], fun.d3FdSigma3->ababbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdSigma3->abbbbb[i], fun.d3FdSigma3->abbbbb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdSigma3->bbbbbb[i], fun.d3FdSigma3->bbbbbb[i], 1e-10, 5e-3);
  }
}

/**
 * @test LibXCTest
 * @brief Tests the return values of XCFun vs LibXC
 */
TEST_F(LibXCTest, Comparison_PBE_Restricted_Cartesian) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "sto-6g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 5;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

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
  const DensityOnGrid<RESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 5e-4);
  ASSERT_NEAR(fun.energy, lib.energy, 1e-6);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-5);
    ASSERT_NEAR((*fun.dFdRho)[i], (*lib.dFdRho)[i], 5.0e-5);
    sAssertNear(d[i], d[i], lib.dFdGradRho->x[i], fun.dFdGradRho->x[i], 1e-12, 1e-5);
    sAssertNear(d[i], d[i], lib.dFdGradRho->y[i], fun.dFdGradRho->y[i], 1e-12, 1e-5);
    sAssertNear(d[i], d[i], lib.dFdGradRho->z[i], fun.dFdGradRho->z[i], 1e-12, 1e-5);
    sAssertNear(d[i], d[i], (*lib.d2FdRho2)[i], (*fun.d2FdRho2)[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdRhodGradRho->x[i], fun.d2FdRhodGradRho->x[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdRhodGradRho->y[i], fun.d2FdRhodGradRho->y[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdRhodGradRho->z[i], fun.d2FdRhodGradRho->z[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->xx[i], fun.d2FdGradRho2->xx[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->xy[i], fun.d2FdGradRho2->xy[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->xz[i], fun.d2FdGradRho2->xz[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->yy[i], fun.d2FdGradRho2->yy[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->yz[i], fun.d2FdGradRho2->yz[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], lib.d2FdGradRho2->zz[i], fun.d2FdGradRho2->zz[i], 1e-10, 5e-4);
    sAssertNear(d[i], d[i], (*lib.d3FdRho3)[i], (*fun.d3FdRho3)[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRho2dGradRho->x[i], fun.d3FdRho2dGradRho->x[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRho2dGradRho->y[i], fun.d3FdRho2dGradRho->y[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRho2dGradRho->z[i], fun.d3FdRho2dGradRho->z[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->xx[i], fun.d3FdRhodGradRho2->xx[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->xy[i], fun.d3FdRhodGradRho2->xy[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->xz[i], fun.d3FdRhodGradRho2->xz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->yy[i], fun.d3FdRhodGradRho2->yy[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->yz[i], fun.d3FdRhodGradRho2->yz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdRhodGradRho2->zz[i], fun.d3FdRhodGradRho2->zz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xxx[i], fun.d3FdGradRho3->xxx[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xxy[i], fun.d3FdGradRho3->xxy[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xxz[i], fun.d3FdGradRho3->xxz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xyy[i], fun.d3FdGradRho3->xyy[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xyz[i], fun.d3FdGradRho3->xyz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->xzz[i], fun.d3FdGradRho3->xzz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->yyy[i], fun.d3FdGradRho3->yyy[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->yyz[i], fun.d3FdGradRho3->yyz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->yzz[i], fun.d3FdGradRho3->yzz[i], 1e-10, 5e-3);
    sAssertNear(d[i], d[i], lib.d3FdGradRho3->zzz[i], fun.d3FdGradRho3->zzz[i], 1e-10, 5e-3);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
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
  const DensityOnGrid<UNRESTRICTED>& d = densOnGridController->getDensityOnGrid();
  auto functional = CompositeFunctionals::resolveFunctional(CompositeFunctionals::FUNCTIONALS::PBE);

  // Calculate data
  // (Second derivatives with libxc for this test included 'nan' values, they have been disabled)
  auto lib = libxc.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);
  auto fun = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, functional, densOnGridController, 3);

  // Compare
  ASSERT_NEAR(fun.epuv->sum(), lib.epuv->sum(), 0.2);
  ASSERT_NEAR(fun.energy, lib.energy, 1e-5);
  for (unsigned int i = 0; i < densOnGridController->getNGridPoints(); i++) {
    // Check all results using absolute values
    ASSERT_NEAR((*fun.epuv)[i], (*lib.epuv)[i], 1.0e-3);
    ASSERT_NEAR(fun.dFdRho->alpha[i], lib.dFdRho->alpha[i], 1.0e-4);
    ASSERT_NEAR(fun.dFdRho->beta[i], lib.dFdRho->beta[i], 1.0e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.dFdGradRho->x.alpha[i], fun.dFdGradRho->x.alpha[i], 1e-12, 1e-5);
    sAssertNear(d.beta[i], d.beta[i], lib.dFdGradRho->x.beta[i], fun.dFdGradRho->x.beta[i], 1e-12, 1e-5);
    sAssertNear(d.alpha[i], d.alpha[i], lib.dFdGradRho->y.alpha[i], fun.dFdGradRho->y.alpha[i], 1e-12, 1e-5);
    sAssertNear(d.beta[i], d.beta[i], lib.dFdGradRho->y.beta[i], fun.dFdGradRho->y.beta[i], 1e-12, 1e-5);
    sAssertNear(d.alpha[i], d.alpha[i], lib.dFdGradRho->z.alpha[i], fun.dFdGradRho->z.alpha[i], 1e-12, 1e-5);
    sAssertNear(d.beta[i], d.beta[i], lib.dFdGradRho->z.beta[i], fun.dFdGradRho->z.beta[i], 1e-12, 1e-5);

    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRho2->aa[i], fun.d2FdRho2->aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRho2->ab[i], fun.d2FdRho2->ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRho2->bb[i], fun.d2FdRho2->bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRhodGradRho->x.aa[i], fun.d2FdRhodGradRho->x.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodGradRho->x.ab[i], fun.d2FdRhodGradRho->x.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRhodGradRho->x.bb[i], fun.d2FdRhodGradRho->x.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRhodGradRho->y.aa[i], fun.d2FdRhodGradRho->y.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodGradRho->y.ab[i], fun.d2FdRhodGradRho->y.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRhodGradRho->y.bb[i], fun.d2FdRhodGradRho->y.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdRhodGradRho->z.aa[i], fun.d2FdRhodGradRho->z.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdRhodGradRho->z.ab[i], fun.d2FdRhodGradRho->z.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdRhodGradRho->z.bb[i], fun.d2FdRhodGradRho->z.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->xx.aa[i], fun.d2FdGradRho2->xx.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->xx.ab[i], fun.d2FdGradRho2->xx.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->xx.bb[i], fun.d2FdGradRho2->xx.bb[i], 1e-10, 5e-4);

    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->xy.aa[i], fun.d2FdGradRho2->xy.aa[i], 1e-10, 1e-5);

    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->xy.ab[i], fun.d2FdGradRho2->xy.ab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->xy.ba[i], fun.d2FdGradRho2->xy.ba[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->xy.bb[i], fun.d2FdGradRho2->xy.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->xz.aa[i], fun.d2FdGradRho2->xz.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->xz.ab[i], fun.d2FdGradRho2->xz.ab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->xz.ba[i], fun.d2FdGradRho2->xz.ba[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->xz.bb[i], fun.d2FdGradRho2->xz.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->yy.aa[i], fun.d2FdGradRho2->yy.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->yy.ab[i], fun.d2FdGradRho2->yy.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->yy.bb[i], fun.d2FdGradRho2->yy.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->yz.aa[i], fun.d2FdGradRho2->yz.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->yz.ab[i], fun.d2FdGradRho2->yz.ab[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->yz.ba[i], fun.d2FdGradRho2->yz.ba[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->yz.bb[i], fun.d2FdGradRho2->yz.bb[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d2FdGradRho2->zz.aa[i], fun.d2FdGradRho2->zz.aa[i], 1e-10, 5e-4);
    sAssertNear(d.alpha[i], d.beta[i], lib.d2FdGradRho2->zz.ab[i], fun.d2FdGradRho2->zz.ab[i], 1e-10, 5e-4);
    sAssertNear(d.beta[i], d.beta[i], lib.d2FdGradRho2->zz.bb[i], fun.d2FdGradRho2->zz.bb[i], 1e-10, 5e-4);

    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho3->aaa[i], fun.d3FdRho3->aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->aab[i], fun.d3FdRho3->aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho3->abb[i], fun.d3FdRho3->abb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho3->bbb[i], fun.d3FdRho3->bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho2dGradRho->x.aaa[i], fun.d3FdRho2dGradRho->x.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->x.aab[i], fun.d3FdRho2dGradRho->x.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->x.aba[i], fun.d3FdRho2dGradRho->x.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->x.abb[i], fun.d3FdRho2dGradRho->x.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->x.bba[i], fun.d3FdRho2dGradRho->x.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho2dGradRho->x.bbb[i], fun.d3FdRho2dGradRho->x.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho2dGradRho->y.aaa[i], fun.d3FdRho2dGradRho->y.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->y.aab[i], fun.d3FdRho2dGradRho->y.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->y.aba[i], fun.d3FdRho2dGradRho->y.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->y.abb[i], fun.d3FdRho2dGradRho->y.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->y.bba[i], fun.d3FdRho2dGradRho->y.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho2dGradRho->y.bbb[i], fun.d3FdRho2dGradRho->y.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRho2dGradRho->z.aaa[i], fun.d3FdRho2dGradRho->z.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->z.aab[i], fun.d3FdRho2dGradRho->z.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->z.aba[i], fun.d3FdRho2dGradRho->z.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->z.abb[i], fun.d3FdRho2dGradRho->z.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRho2dGradRho->z.bba[i], fun.d3FdRho2dGradRho->z.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRho2dGradRho->z.bbb[i], fun.d3FdRho2dGradRho->z.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->xx.aaa[i], fun.d3FdRhodGradRho2->xx.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xx.aab[i], fun.d3FdRhodGradRho2->xx.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xx.abb[i], fun.d3FdRhodGradRho2->xx.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xx.baa[i], fun.d3FdRhodGradRho2->xx.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xx.bab[i], fun.d3FdRhodGradRho2->xx.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->xx.bbb[i], fun.d3FdRhodGradRho2->xx.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->xy.aaa[i], fun.d3FdRhodGradRho2->xy.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.aab[i], fun.d3FdRhodGradRho2->xy.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.aba[i], fun.d3FdRhodGradRho2->xy.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.abb[i], fun.d3FdRhodGradRho2->xy.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.baa[i], fun.d3FdRhodGradRho2->xy.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.bab[i], fun.d3FdRhodGradRho2->xy.bab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xy.bba[i], fun.d3FdRhodGradRho2->xy.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->xy.bbb[i], fun.d3FdRhodGradRho2->xy.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->xz.aaa[i], fun.d3FdRhodGradRho2->xz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.aab[i], fun.d3FdRhodGradRho2->xz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.aba[i], fun.d3FdRhodGradRho2->xz.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.abb[i], fun.d3FdRhodGradRho2->xz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.baa[i], fun.d3FdRhodGradRho2->xz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.bab[i], fun.d3FdRhodGradRho2->xz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->xz.bba[i], fun.d3FdRhodGradRho2->xz.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->xz.bbb[i], fun.d3FdRhodGradRho2->xz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->yy.aaa[i], fun.d3FdRhodGradRho2->yy.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yy.aab[i], fun.d3FdRhodGradRho2->yy.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yy.abb[i], fun.d3FdRhodGradRho2->yy.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yy.baa[i], fun.d3FdRhodGradRho2->yy.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yy.bab[i], fun.d3FdRhodGradRho2->yy.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->yy.bbb[i], fun.d3FdRhodGradRho2->yy.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->yz.aaa[i], fun.d3FdRhodGradRho2->yz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.aab[i], fun.d3FdRhodGradRho2->yz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.aba[i], fun.d3FdRhodGradRho2->yz.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.abb[i], fun.d3FdRhodGradRho2->yz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.baa[i], fun.d3FdRhodGradRho2->yz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.bab[i], fun.d3FdRhodGradRho2->yz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->yz.bba[i], fun.d3FdRhodGradRho2->yz.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->yz.bbb[i], fun.d3FdRhodGradRho2->yz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdRhodGradRho2->zz.aaa[i], fun.d3FdRhodGradRho2->zz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->zz.aab[i], fun.d3FdRhodGradRho2->zz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->zz.abb[i], fun.d3FdRhodGradRho2->zz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->zz.baa[i], fun.d3FdRhodGradRho2->zz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdRhodGradRho2->zz.bab[i], fun.d3FdRhodGradRho2->zz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdRhodGradRho2->zz.bbb[i], fun.d3FdRhodGradRho2->zz.bbb[i], 1e-10, 5e-3);

    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xxx.aaa[i], fun.d3FdGradRho3->xxx.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxx.aab[i], fun.d3FdGradRho3->xxx.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxx.abb[i], fun.d3FdGradRho3->xxx.abb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xxx.bbb[i], fun.d3FdGradRho3->xxx.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xxy.aaa[i], fun.d3FdGradRho3->xxy.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxy.aab[i], fun.d3FdGradRho3->xxy.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxy.aba[i], fun.d3FdGradRho3->xxy.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxy.abb[i], fun.d3FdGradRho3->xxy.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxy.bba[i], fun.d3FdGradRho3->xxy.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xxy.bbb[i], fun.d3FdGradRho3->xxy.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xxz.aaa[i], fun.d3FdGradRho3->xxz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxz.aab[i], fun.d3FdGradRho3->xxz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxz.aba[i], fun.d3FdGradRho3->xxz.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxz.abb[i], fun.d3FdGradRho3->xxz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xxz.bba[i], fun.d3FdGradRho3->xxz.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xxz.bbb[i], fun.d3FdGradRho3->xxz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xyy.aaa[i], fun.d3FdGradRho3->xyy.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyy.aab[i], fun.d3FdGradRho3->xyy.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyy.abb[i], fun.d3FdGradRho3->xyy.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyy.baa[i], fun.d3FdGradRho3->xyy.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyy.bab[i], fun.d3FdGradRho3->xyy.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xyy.bbb[i], fun.d3FdGradRho3->xyy.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xyz.aaa[i], fun.d3FdGradRho3->xyz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.aab[i], fun.d3FdGradRho3->xyz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.aba[i], fun.d3FdGradRho3->xyz.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.abb[i], fun.d3FdGradRho3->xyz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.baa[i], fun.d3FdGradRho3->xyz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.bab[i], fun.d3FdGradRho3->xyz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xyz.bba[i], fun.d3FdGradRho3->xyz.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xyz.bbb[i], fun.d3FdGradRho3->xyz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->xzz.aaa[i], fun.d3FdGradRho3->xzz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xzz.aab[i], fun.d3FdGradRho3->xzz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xzz.abb[i], fun.d3FdGradRho3->xzz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xzz.baa[i], fun.d3FdGradRho3->xzz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->xzz.bab[i], fun.d3FdGradRho3->xzz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->xzz.bbb[i], fun.d3FdGradRho3->xzz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->yyy.aaa[i], fun.d3FdGradRho3->yyy.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyy.aab[i], fun.d3FdGradRho3->yyy.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyy.abb[i], fun.d3FdGradRho3->yyy.abb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->yyy.bbb[i], fun.d3FdGradRho3->yyy.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->yyz.aaa[i], fun.d3FdGradRho3->yyz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyz.aab[i], fun.d3FdGradRho3->yyz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyz.aba[i], fun.d3FdGradRho3->yyz.aba[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyz.abb[i], fun.d3FdGradRho3->yyz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yyz.bba[i], fun.d3FdGradRho3->yyz.bba[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->yyz.bbb[i], fun.d3FdGradRho3->yyz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->yzz.aaa[i], fun.d3FdGradRho3->yzz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yzz.aab[i], fun.d3FdGradRho3->yzz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yzz.abb[i], fun.d3FdGradRho3->yzz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yzz.baa[i], fun.d3FdGradRho3->yzz.baa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->yzz.bab[i], fun.d3FdGradRho3->yzz.bab[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->yzz.bbb[i], fun.d3FdGradRho3->yzz.bbb[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.alpha[i], lib.d3FdGradRho3->zzz.aaa[i], fun.d3FdGradRho3->zzz.aaa[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->zzz.aab[i], fun.d3FdGradRho3->zzz.aab[i], 1e-10, 5e-3);
    sAssertNear(d.alpha[i], d.beta[i], lib.d3FdGradRho3->zzz.abb[i], fun.d3FdGradRho3->zzz.abb[i], 1e-10, 5e-3);
    sAssertNear(d.beta[i], d.beta[i], lib.d3FdGradRho3->zzz.bbb[i], fun.d3FdGradRho3->zzz.bbb[i], 1e-10, 5e-3);
  }
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
