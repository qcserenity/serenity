/**
 * @file SAOPPotential_test.cpp
 *
 * @date   Aug 4, 2017
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
#include "potentials/SAOPPotential.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class SAOPPotentialTest : public ::testing::Test {
 protected:
  SAOPPotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Ne2_6_31Gs)) {
  }
  virtual ~SAOPPotentialTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// the system
  std::shared_ptr<SystemController> systemController;
};

/**
 * @test SAOPPotentialTest
 * @brief Test the SAOP-potential of a neon dimer (restricted).
 */
TEST_F(SAOPPotentialTest, Ne_dimer_rFockMatrix) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  auto basis = systemController->getBasisController();
  auto grid = systemController->getGridController();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(systemController->getSettings(), basis, grid);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densMatrixCon = systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  auto densOnGridCon = std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalculator, densMatrixCon);
  SAOPPotential<SCFMode> saopCalc(basis, systemController, densOnGridCon);
  auto F = saopCalc.getMatrix();
  /*
   * Warning! This result was obtained with Serenity with a working SAOP version (18.09.2017)
   * and not further verified.
   */
  EXPECT_NEAR(-32.846378197834504, F.trace(), 1e-5);
  EXPECT_NEAR(-25.079698138378454, saopCalc.getEnergy(densMatrixCon->getDensityMatrix()), 1e-5);
}

/**
 * @test SAOPPotentialTest
 * @brief Test the SAOP-potential of a neon dimer (unrestricted).
 */
TEST_F(SAOPPotentialTest, Ne_dimer_uFockMatrix) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;
  auto basis = systemController->getBasisController();
  auto grid = systemController->getGridController();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(systemController->getSettings(), basis, grid);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densMatrixCon = systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  auto densOnGridCon = std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densOnGridCalculator, densMatrixCon);
  SAOPPotential<SCFMode> saopCalc(basis, systemController, densOnGridCon);
  auto F(std::move(saopCalc.getMatrix()));
  /*
   * Warning! This result was obtained with Serenity with a working SAOP version (18.09.2017)
   * and not further verified.
   */
  EXPECT_NEAR(-32.846378274630261, F.alpha.trace(), 1e-5);
  EXPECT_NEAR(-32.846378274630261, F.beta.trace(), 1e-5);
  EXPECT_NEAR(-25.079698138378454, saopCalc.getEnergy(densMatrixCon->getDensityMatrix()), 1e-5);
}

} /* namespace Serenity */
