/**
 * @file   DensityAdder_test.cpp
 * @author M. Boeckers
 *
 * @date   06. April 2018
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
#include "data/grid/DensityAdder.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class DensityAdderTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(DensityAdderTest, RDensityAdder) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto D = sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(sys->getSettings(), sys->getBasisController(), sys->getGridController());
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(
      basisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
  auto densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>(
      sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>(),
      sys->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>());
  auto densityMatrixDensityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>(densityOnGridCalculator,
                                                                                             densityMatrixController);
  auto p = densityMatrixDensityOnGridController->getDensityOnGrid();
  DensityOnGrid<Options::SCF_MODES::RESTRICTED> ptest(sys->getGridController());
  DensityAdder<Options::SCF_MODES::RESTRICTED>::add(ptest, D, sys->getBasisController(), 0.0);
  EXPECT_NEAR((p - ptest).norm(), 0.0, 1.0e-6);
}

TEST_F(DensityAdderTest, UDensityAdder) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto D = sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(sys->getSettings(), sys->getBasisController(), sys->getGridController());
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
  auto densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
      sys->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>(),
      sys->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
  auto densityMatrixDensityOnGridController =
      std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(densityOnGridCalculator,
                                                                                               densityMatrixController);
  auto& p = densityMatrixDensityOnGridController->getDensityOnGrid();
  DensityOnGrid<Options::SCF_MODES::UNRESTRICTED> ptest(sys->getGridController());
  DensityAdder<Options::SCF_MODES::UNRESTRICTED>::add(ptest, D, sys->getBasisController(), 0.0);
  EXPECT_NEAR((p.alpha - ptest.alpha).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((p.beta - ptest.beta).norm(), 0.0, 1.0e-6);
}

} // namespace Serenity
