/**
 * @file DensityMatrixDensityOnGridController_test.cpp
 *
 * @date 26 Mar 2020
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
#include "data/grid/DensityMatrixDensityOnGridController.h" //To be tested
#include "data/ElectronicStructure.h"                       //getDensityMatrixController.
#include "data/grid/DensityOnGridFactory.h"                 //DensityMatrixDensityOnGridController construction.
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test supply.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
namespace Serenity {

class DensityMatrixDensityOnGridControllerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(DensityMatrixDensityOnGridControllerTest, setHighestDerivativeViaGetHessian) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto densityOnGridController =
      DensityOnGridFactory<RESTRICTED>::produce(system->getElectronicStructure<RESTRICTED>()->getDensityMatrixController(),
                                                system->getGridController(), 1, system->getSettings());
  auto densityHessian = densityOnGridController->getDensityHessianOnGrid();
  auto refDensity = densityOnGridController->getDensityOnGrid();

  auto refDensityOnGridController =
      DensityOnGridFactory<RESTRICTED>::produce(system->getElectronicStructure<RESTRICTED>()->getDensityMatrixController(),
                                                system->getGridController(), 2, system->getSettings());

  auto refHessian = refDensityOnGridController->getDensityHessianOnGrid();
  double diff =
      (refHessian.xx - densityHessian.xx).array().abs().sum() + (refHessian.xy - densityHessian.xy).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-9);
  densityOnGridController->setHighestDerivative(0);
  auto density = densityOnGridController->getDensityOnGrid();
  diff = (density - refDensity).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-9);
}

} /* namespace Serenity */
