/**
 * @file OptEffPotential_test.cpp
 *
 * @date Dec 19, 2016
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
#include "potentials/OptEffPotential.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/GridPotential.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/FockMatrix.h"
#include "geometry/Geometry.h"
#include "integrals/OneIntControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class OptEffPotentialTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    auto& libint = Libint::getInstance();
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
    auto& libint = Libint::getInstance();
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }
};

TEST_F(OptEffPotentialTest, oep) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);

  auto gridController = system->getGridController();
  auto basisController = system->getBasisController();
  auto basisFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  auto oneEIntController = system->getOneElectronIntegralController();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(
      basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  auto densMat = system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
  const auto& targetDens = densOnGridCalculator->calcDensityOnGrid(densMat);
  auto nOccOrbs = system->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto resultOrbitals = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(basisController, 0);

  GridPotential<Options::SCF_MODES::RESTRICTED> initialGuess(gridController);
  CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::RESTRICTED>(initialGuess, densMat);
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(initialGuess, system->getGeometry()->getAtoms());

  FockMatrix<Options::SCF_MODES::RESTRICTED> initialGuessMat(basisController);

  ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED> scalarOpToMat(basisFunctionOnGridController, 0.0);
  scalarOpToMat.addScalarOperatorToMatrix(initialGuessMat, initialGuess);

  OptEffPotential<Options::SCF_MODES::RESTRICTED> oepCalc(basisFuncOnGridController, basisFuncOnGridController,
                                                          oneEIntController, densOnGridCalculator, nOccOrbs, 1e-5, 1e-5, 0.2);

  oepCalc.calculateOEP(targetDens, resultOrbitals, initialGuess);

  initialGuessMat += oepCalc.getMatrix();

  oepCalc.calculateOEPCarter(densMat, resultOrbitals, initialGuessMat);

  auto& oep = oepCalc.getPotentialOnGrid();

  DensityMatrixController<Options::SCF_MODES::RESTRICTED> recDMatController(
      resultOrbitals,
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()->getOccupations());

  const auto& recDens = densOnGridCalculator->calcDensityOnGrid(recDMatController.getDensityMatrix());

  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd> dummy(
      Eigen::VectorXd::Zero(basisController->getNBasisFunctions()));

  auto grad = oepCalc.getGradient(targetDens, recDens, dummy);

  /*
   * Gradients should be nearly zero
   */

  initialGuess += oep;

  auto densDiff = oepCalc.calculateOEPLB(targetDens, resultOrbitals, initialGuess);

  /*
   * Density difference should be low
   */

  EXPECT_NEAR(0.0, densDiff, 5.0e-3);
}

TEST_F(OptEffPotentialTest, oepUnres) {
  // get test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);

  auto gridController = system->getGridController();
  auto basisController = system->getBasisController();
  auto basisFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  auto oneEIntController = system->getOneElectronIntegralController();
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(system->getSettings(), basisController, gridController);
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, system->getSettings().grid.blockAveThreshold);
  auto densMat = system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix();
  const auto& targetDens = densOnGridCalculator->calcDensityOnGrid(densMat);
  auto nOccOrbs = system->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto resultOrbitals = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(basisController, 0);

  GridPotential<Options::SCF_MODES::UNRESTRICTED> initialGuess(gridController);
  GridPotential<Options::SCF_MODES::RESTRICTED> tmp(gridController);
  CoulombPotentialOnGridCalculator::calculateElectronElectron<Options::SCF_MODES::UNRESTRICTED>(tmp, densMat);
  CoulombPotentialOnGridCalculator::calculateElectronNuclei(tmp, system->getGeometry()->getAtoms());

  for_spin(initialGuess) {
    initialGuess_spin += tmp;
  };

  OptEffPotential<Options::SCF_MODES::UNRESTRICTED> oepCalc(basisFuncOnGridController, basisFuncOnGridController,
                                                            oneEIntController, densOnGridCalculator, nOccOrbs, 1e-5,
                                                            1e-5, 0.2);

  oepCalc.calculateOEP(targetDens, resultOrbitals, initialGuess);

  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> recDMatController(
      resultOrbitals,
      system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController()->getOccupations());

  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> dummy(
      Eigen::VectorXd::Zero(basisController->getNBasisFunctions()));

  /*
   * Gradients should be nearly zero, for unrestricted the threshold should be less strong
   */

  auto& oep = oepCalc.getPotentialOnGrid();

  initialGuess += oep;

  const auto& densDiff = oepCalc.calculateOEPLB(targetDens, resultOrbitals, initialGuess);

  /*
   * Density difference should be low
   */

  EXPECT_NEAR(0.0, densDiff.alpha, 5.0e-3);
  EXPECT_NEAR(0.0, densDiff.beta, 5.0e-3);
}

} /*namespace Serenity*/
