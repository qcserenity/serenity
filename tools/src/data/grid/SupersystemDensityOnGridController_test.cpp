/**
 * @file SupersystemDensityOnGridController_test.cpp
 *
 * @date Oct 12, 2017
 * @author Jan Unsleber
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
#include "data/grid/SupersystemDensityOnGridController.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "geometry/AtomTypeFactory.h"
#include "grid/GridController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class SupersystemDensityOnGridControllerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests SupersystemDensityOnGridController.h/.cpp.
 */
TEST_F(SupersystemDensityOnGridControllerTest, Restricted) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  auto dMat = system->getElectronicStructure<RESTRICTED>()->getDensityMatrixController();

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      system->getSettings(), dMat->getDensityMatrix().getBasisController(), system->getGridController());
  basisFunctionOnGridController->setHighestDerivative(2);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  auto subsys = std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(densityOnGridCalculator, dMat, 2);
  std::vector<std::shared_ptr<DensityOnGridController<RESTRICTED>>> vec;
  vec.push_back(subsys);
  vec.push_back(subsys);
  SupersystemDensityOnGridController<RESTRICTED> super(vec);
  const auto& densSuper = super.getDensityOnGrid();
  const auto& densSub = subsys->getDensityOnGrid();
  unsigned int nPoints = densSuper.size();
  // density
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(densSuper[i], densSub[i] + densSub[i]);
  }
  // gradient
  const auto& gradSuper = super.getDensityGradientOnGrid();
  const auto& gradSub = subsys->getDensityGradientOnGrid();
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(gradSuper.x[i], gradSub.x[i] + gradSub.x[i]);
    EXPECT_EQ(gradSuper.y[i], gradSub.y[i] + gradSub.y[i]);
    EXPECT_EQ(gradSuper.z[i], gradSub.z[i] + gradSub.z[i]);
  }
  // hessian
  const auto& hessSuper = super.getDensityHessianOnGrid();
  const auto& hessSub = subsys->getDensityHessianOnGrid();
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(hessSuper.xx[i], hessSub.xx[i] + hessSub.xx[i]);
    EXPECT_EQ(hessSuper.xy[i], hessSub.xy[i] + hessSub.xy[i]);
    EXPECT_EQ(hessSuper.xz[i], hessSub.xz[i] + hessSub.xz[i]);
    EXPECT_EQ(hessSuper.yy[i], hessSub.yy[i] + hessSub.yy[i]);
    EXPECT_EQ(hessSuper.yz[i], hessSub.yz[i] + hessSub.yz[i]);
    EXPECT_EQ(hessSuper.zz[i], hessSub.zz[i] + hessSub.zz[i]);
  }
}

/**
 * @test
 * @brief Tests SupersystemDensityOnGridController.h/.cpp.
 */
TEST_F(SupersystemDensityOnGridControllerTest, Unrestricted) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  auto dMat = system->getElectronicStructure<UNRESTRICTED>()->getDensityMatrixController();

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      system->getSettings(), dMat->getDensityMatrix().getBasisController(), system->getGridController());
  basisFunctionOnGridController->setHighestDerivative(2);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(
      basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  auto subsys = std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(densityOnGridCalculator, dMat, 2);
  std::vector<std::shared_ptr<DensityOnGridController<UNRESTRICTED>>> vec;
  vec.push_back(subsys);
  vec.push_back(subsys);
  SupersystemDensityOnGridController<UNRESTRICTED> super(vec);
  const auto& densSuper = super.getDensityOnGrid();
  const auto& densSub = subsys->getDensityOnGrid();
  unsigned int nPoints = densSuper.alpha.size();
  // density
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(densSuper.alpha[i], densSub.alpha[i] + densSub.alpha[i]);
    EXPECT_EQ(densSuper.beta[i], densSub.beta[i] + densSub.beta[i]);
  }
  // gradient
  const auto& gradSuper = super.getDensityGradientOnGrid();
  const auto& gradSub = subsys->getDensityGradientOnGrid();
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(gradSuper.x.alpha[i], gradSub.x.alpha[i] + gradSub.x.alpha[i]);
    EXPECT_EQ(gradSuper.y.alpha[i], gradSub.y.alpha[i] + gradSub.y.alpha[i]);
    EXPECT_EQ(gradSuper.z.alpha[i], gradSub.z.alpha[i] + gradSub.z.alpha[i]);
    EXPECT_EQ(gradSuper.x.beta[i], gradSub.x.beta[i] + gradSub.x.beta[i]);
    EXPECT_EQ(gradSuper.y.beta[i], gradSub.y.beta[i] + gradSub.y.beta[i]);
    EXPECT_EQ(gradSuper.z.beta[i], gradSub.z.beta[i] + gradSub.z.beta[i]);
  }
  // hessian
  const auto& hessSuper = super.getDensityHessianOnGrid();
  const auto& hessSub = subsys->getDensityHessianOnGrid();
  for (unsigned int i = 0; i < nPoints; i++) {
    EXPECT_EQ(hessSuper.xx.alpha[i], hessSub.xx.alpha[i] + hessSub.xx.alpha[i]);
    EXPECT_EQ(hessSuper.xy.alpha[i], hessSub.xy.alpha[i] + hessSub.xy.alpha[i]);
    EXPECT_EQ(hessSuper.xz.alpha[i], hessSub.xz.alpha[i] + hessSub.xz.alpha[i]);
    EXPECT_EQ(hessSuper.yy.alpha[i], hessSub.yy.alpha[i] + hessSub.yy.alpha[i]);
    EXPECT_EQ(hessSuper.yz.alpha[i], hessSub.yz.alpha[i] + hessSub.yz.alpha[i]);
    EXPECT_EQ(hessSuper.zz.alpha[i], hessSub.zz.alpha[i] + hessSub.zz.alpha[i]);
    EXPECT_EQ(hessSuper.xx.beta[i], hessSub.xx.beta[i] + hessSub.xx.beta[i]);
    EXPECT_EQ(hessSuper.xy.beta[i], hessSub.xy.beta[i] + hessSub.xy.beta[i]);
    EXPECT_EQ(hessSuper.xz.beta[i], hessSub.xz.beta[i] + hessSub.xz.beta[i]);
    EXPECT_EQ(hessSuper.yy.beta[i], hessSub.yy.beta[i] + hessSub.yy.beta[i]);
    EXPECT_EQ(hessSuper.yz.beta[i], hessSub.yz.beta[i] + hessSub.yz.beta[i]);
    EXPECT_EQ(hessSuper.zz.beta[i], hessSub.zz.beta[i] + hessSub.zz.beta[i]);
  }
}

} /* namespace Serenity */
