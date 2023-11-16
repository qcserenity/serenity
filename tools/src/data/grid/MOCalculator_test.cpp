/**
 * @file MOCalculator_test.cpp
 *
 * @date Sep 13, 2016
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
#include "data/grid/MOCalculator.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/matrices/CoefficientMatrix.h"
#include "grid/GridController.h"
#include "math/Matrix.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class MOCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(MOCalculatorTest, MOValidate) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);

  // get some useful data
  auto basisController = system->getBasisController();
  auto gridController = system->getGridController();
  auto& weights = gridController->getWeights();
  auto nGridPts = gridController->getNGridPoints();
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  /*
   * make sure to get the orbitals via the ElectronicStructure so that an
   * electronic structure calculation is performed if no orbitals are
   * available
   */
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();

  // the MOCalculator to test
  MOCalculator moCalc(basFuncOnGridController);

  /*
   * Get the MOs on a grid. Call calcAllMOValues since it calls the other two
   * functions and tests them as well.
   */
  auto allMos = moCalc.calcAllMOValuesOnGrid<Options::SCF_MODES::RESTRICTED>(coeffs, 1e-9);

  // Build integrals <1|1> and <1|2>
  double intMo = 0.0;
  double intMoMo = 0.0;
  for (unsigned int gridPt = 0; gridPt < nGridPts; gridPt++) {
    intMo += (weights[gridPt] * allMos(gridPt, 0) * allMos(gridPt, 0));
    intMoMo += (weights[gridPt] * (allMos(gridPt, 0) * allMos(gridPt, 1)));
  }

  // Orbitals should be orthonormal
  EXPECT_NEAR(1.0, intMo, 1.0e-6);
  EXPECT_NEAR(0.0, intMoMo, 1.0e-6);
}

TEST_F(MOCalculatorTest, MOValidate_Unrestricted) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);

  // get some useful data
  auto basisController = system->getBasisController();
  auto gridController = system->getGridController();
  auto& weights = gridController->getWeights();
  auto nGridPts = gridController->getNGridPoints();
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  /*
   * make sure to get the orbitals via the ElectronicStructure so that an
   * electronic structure calculation is performed if no orbitals are
   * available
   */
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coeffs =
      system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();

  // the MOCalculator to test
  MOCalculator moCalc(basFuncOnGridController);

  /*
   * Get the MOs on a grid. Call calcAllMOValues since it calls the other two
   * functions and tests them as well.
   */
  auto allMos = moCalc.calcAllMOValuesOnGrid<Options::SCF_MODES::UNRESTRICTED>(coeffs);

  // Build integrals <1|1> and <1|2>
  double intMo = 0.0;
  double intMoMo = 0.0;
  for (unsigned int gridPt = 0; gridPt < nGridPts; gridPt++) {
    intMo += (weights[gridPt] * allMos.alpha(gridPt, 0) * allMos.alpha(gridPt, 0));
    intMoMo += (weights[gridPt] * (allMos.beta(gridPt, 0) * allMos.beta(gridPt, 1)));
  }

  // Orbitals should be orthonormal
  EXPECT_NEAR(1.0, intMo, 1.0e-6);
  EXPECT_NEAR(0.0, intMoMo, 1.0e-6);

  // TODO Test kinetic energy density... but how?
}

} // namespace Serenity
