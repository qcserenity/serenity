/**
 * @file   CM5PopulationCalculator_test.cpp
 *
 * @date   Oct 7, 2024
 * @author Thorben Wiegmann
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
#include "analysis/populationAnalysis/CM5PopulationCalculator.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory>

namespace Serenity {

/**
 * @class CM5PopulationCalculatorTest
 * @brief Sets everything up for the tests of CM5PopulationCalculator.h/.cpp.\n
          All of these results were obtained with Gaussian 16 on October 07, 2024.
 */
class CM5PopulationCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(CM5PopulationCalculatorTest, waterRestricted) {
  /*
   * Set up test environment
   */
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), systemController->getBasisController(), systemController->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(
      densOnGridCalculator, systemController->template getElectronicStructure<RESTRICTED>()->getDensityMatrixController());
  std::shared_ptr<HirshfeldPopulationCalculator<RESTRICTED>> hirshFeld =
      std::make_shared<HirshfeldPopulationCalculator<RESTRICTED>>(systemController, densOnGridController);
  CM5PopulationCalculator<RESTRICTED> calculator(systemController, hirshFeld);
  auto populations = calculator.getAtomPopulations();
  /*
  /*
   * Check results
   */
  EXPECT_NEAR(1 - populations[0], 0.331197, 1e-2);
  EXPECT_NEAR(1 - populations[1], 0.331197, 1e-2);
  EXPECT_NEAR(8 - populations[2], -0.662394, 1e-2);
  // the HirshfeldPopulationCalculator performs SCF_INPLACE for the individual atoms -> remove their system folders
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "O_FREE/", "O_FREE");
}

TEST_F(CM5PopulationCalculatorTest, waterUnrestricted) {
  /*
   * Set up test environment
   */
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), systemController->getBasisController(), systemController->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(
      densOnGridCalculator, systemController->template getElectronicStructure<UNRESTRICTED>()->getDensityMatrixController());
  std::shared_ptr<HirshfeldPopulationCalculator<UNRESTRICTED>> hirshFeld =
      std::make_shared<HirshfeldPopulationCalculator<UNRESTRICTED>>(systemController, densOnGridController);
  CM5PopulationCalculator<UNRESTRICTED> calculator(systemController, hirshFeld);
  auto populations = calculator.getAtomPopulations();
  const Eigen::VectorXd totPops = populations.total();
  /*
  /*
   * Check results
   */
  EXPECT_NEAR(1 - totPops[0], 0.331197, 1e-2);
  EXPECT_NEAR(1 - totPops[1], 0.331197, 1e-2);
  EXPECT_NEAR(8 - totPops[2], -0.662394, 1e-2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "O_FREE/", "O_FREE");
}

} // namespace Serenity