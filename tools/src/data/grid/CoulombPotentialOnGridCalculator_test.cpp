/**
 * @file CoulombPotentialOnGridCalculator_test.cpp
 *
 * @date Mar 31, 2016
 * @author David Schnieders
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
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "data/ElectronicStructure.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridFactory.h"
#include "geometry/Geometry.h"
#include "grid/GridController.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class CoulombPotentialOnGridCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(CoulombPotentialOnGridCalculatorTest, completeCoul_Restricted) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  // get some useful data
  auto atoms = system->getGeometry()->getAtoms();
  auto basisController = system->getBasisController();
  auto gridController = system->getGridController();
  GridPotential<RESTRICTED> coulGrid(gridController);
  /*
   * make sure to get the orbitals via the ElectronicStructure so that an
   * electronic structure calculation is performed if no orbitals are
   * available
   */
  auto densMat = system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();

  // the CoulombPotentialOnGridCalculator to test
  CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::RESTRICTED>(coulGrid, densMat);
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(coulGrid, atoms);

  // Orbitals should be orthonormal
  EXPECT_NEAR(0.00024478514611243896, coulGrid[22], 1.0e-6);
  EXPECT_NEAR(2.0298323251469186e-05, coulGrid[1], 1.0e-6);
  EXPECT_NEAR(5.9034570994223601e-05, coulGrid[10], 1.0e-6);
  EXPECT_NEAR(0.00018177569872274146, coulGrid[100], 1.0e-6);
}

TEST_F(CoulombPotentialOnGridCalculatorTest, completeCoul_Unrestricted) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  // get some useful data
  auto atoms = system->getGeometry()->getAtoms();
  auto basisController = system->getBasisController();
  auto gridController = system->getGridController();
  GridPotential<RESTRICTED> coulGrid(gridController);
  /*
   * make sure to get the orbitals via the ElectronicStructure so that an
   * electronic structure calculation is performed if no orbitals are
   * available
   */
  auto densMat = system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix();

  // the CoulombPotentialOnGridCalculator to test
  CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::UNRESTRICTED>(coulGrid, densMat);
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(coulGrid, atoms);

  // Orbitals should be orthonormal
  EXPECT_NEAR(0.00024478514611243896, coulGrid[22], 1.0e-6);
  EXPECT_NEAR(2.0298323251469186e-05, coulGrid[1], 1.0e-6);
  EXPECT_NEAR(5.9034570994223601e-05, coulGrid[10], 1.0e-6);
  EXPECT_NEAR(0.00018177569872274146, coulGrid[100], 1.0e-6);
}

} /* namespace Serenity */
