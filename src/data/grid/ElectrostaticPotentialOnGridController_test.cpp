/**
 * @file ElectrostaticPotentialOnGridController_test.cpp
 *
 * @author Moritz Bensberg
 * @date Aug 6, 2020
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
#include "data/grid/ElectrostaticPotentialOnGridController.h" //To be tested.
#include "data/ElectronicStructure.h"                         //getDensityMatrixController.
#include "geometry/MolecularSurfaceController.h"              //getCavityGridController
#include "system/SystemController.h"                          //Test systems/cavity.
#include "testsupply/SystemController__TEST_SUPPLY.h"         //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class ElectrostaticPotentialOnGridControllerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ElectrostaticPotentialOnGridControllerTest, WaterOnCavity) {
  auto water = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto elecPotOnGridController =
      water->getElectrostaticPotentialOnMolecularSurfaceController<RESTRICTED>(MOLECULAR_SURFACE_TYPES::ACTIVE);

  const GridPotential<RESTRICTED> potential = elecPotOnGridController->getPotential();
  const double testSum = potential.sum();
  EXPECT_NEAR(testSum, -0.80515632923108371, 1e-5);

  elecPotOnGridController->setDiskMode(true);
  const GridPotential<RESTRICTED> potentialFromDisk = elecPotOnGridController->getPotential();
  double diff = (potentialFromDisk - potential).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-6);

  elecPotOnGridController->notify();
  const GridPotential<RESTRICTED> recalculatedPotential = elecPotOnGridController->getPotential();
  diff = (recalculatedPotential - potential).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-6);

  elecPotOnGridController->cleanUpDisk();
  SystemController__TEST_SUPPLY::cleanUp();
}
} /* namespace Serenity */
