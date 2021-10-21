/**
 * @file   HirshfeldPopulationCalculator_test.cpp
 *
 * @date   Jul 4, 2017
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
#include "analysis/populationAnalysis/HirshfeldPopulationCalculator.h"
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

class HirshfeldPopulationCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(HirshfeldPopulationCalculatorTest, populationsRestricted) {
  /*
   * Set up test environment
   */
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), systemController->getBasisController(), systemController->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(
      densOnGridCalculator, systemController->template getElectronicStructure<RESTRICTED>()->getDensityMatrixController());

  HirshfeldPopulationCalculator<RESTRICTED> hirshFeld(systemController, densOnGridController);

  auto populations = hirshFeld.getAtomPopulations();
  /*
   * Check results
   *
   *
   * The promolecular density is generated from spherically averaged atom densities, obtained from
   * scf calculations using the same method as specified for the target system. These scf calculations will
   * be done in-place using quite weak convergence criteria. This may lead to some numerical instabilities
   * in the resulting atom populations (in the 4th decimal). If a higher accuracy is needed, the atom scf
   * convergence criteria need to be tightened (take a look at the AtomDensityGuessCalculator for this).
   */

  EXPECT_NEAR(populations[0], 5.8875, 1e-1);
  EXPECT_NEAR(populations[1], 8.1124, 1e-1);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "C_FREE/", "C_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "O_FREE/", "O_FREE");
}
/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(HirshfeldPopulationCalculatorTest, populationsRestricted2) {
  /*
   * Set up test environment
   */
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), systemController->getBasisController(), systemController->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(
      densOnGridCalculator, systemController->template getElectronicStructure<RESTRICTED>()->getDensityMatrixController());

  HirshfeldPopulationCalculator<RESTRICTED> hirshFeld(systemController, densOnGridController);

  auto populations = hirshFeld.getAtomPopulations();
  /*
   * Check results
   *
   *
   * The promolecular density is generated from spherically averaged atom densities, obtained from
   * scf calculations using the same method as specified for the target system. These scf calculations will
   * be done in-place using quite weak convergence criteria. This may lead to some numerical instabilities
   * in the resulting atom populations (in the 4th decimal). If a higher accuracy is needed, the atom scf
   * convergence criteria need to be tightened (take a look at the AtomDensityGuessCalculator for this).
   */
  EXPECT_NEAR(populations[0], 0.8352, 1e-2);
  EXPECT_NEAR(populations[1], 0.8339, 1e-2);
  EXPECT_NEAR(populations[2], 8.3307, 1e-2);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "O_FREE/", "O_FREE");
}

/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(HirshfeldPopulationCalculatorTest, populationsUnrestricted) {
  /*
   * Set up test environment
   */
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);

  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      systemController->getSettings(), systemController->getBasisController(), systemController->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<UNRESTRICTED>>(
      basisFunctionOnGridController, systemController->getSettings().grid.blockAveThreshold);
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<UNRESTRICTED>>(
      densOnGridCalculator, systemController->template getElectronicStructure<UNRESTRICTED>()->getDensityMatrixController());

  HirshfeldPopulationCalculator<UNRESTRICTED> hirshFeld(systemController, densOnGridController);

  auto populations = hirshFeld.getAtomPopulations();
  const Eigen::VectorXd totalPops = populations.total();
  /*
   * Check results
   *
   *
   * The promolecular density is generated from spherically averaged atom densities, obtained from
   * scf calculations using the same method as specified for the target system. These scf calculations will
   * be done in-place using quite weak convergence criteria. This may lead to some numerical instabilities
   * in the resulting atom populations (in the 4th decimal). If a higher accuracy is needed, the atom scf
   * convergence criteria need to be tightened (take a look at the AtomDensityGuessCalculator for this).
   */
  EXPECT_NEAR(totalPops[0], 0.84174920578360513, 1e-2);
  EXPECT_NEAR(totalPops[1], 0.84163292034736281, 1e-2);
  EXPECT_NEAR(totalPops[2], 8.3166644990895318, 1e-2);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "H_FREE/", "H_FREE");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController->getSystemPath() + "O_FREE/", "O_FREE");
}

} // namespace Serenity
