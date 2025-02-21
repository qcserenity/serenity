/**
 * @file   DensityInitialGuessCalculator.cpp
 * @author Thomas Dresselhaus
 *
 * @date   11. Juli 2014, 17:22
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
/* Include Class Header*/
#include "scf/initialGuess/DensityInitialGuessCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "potentials/ERIPotential.h"
#include "potentials/EffectiveCorePotential.h"
#include "potentials/FuncPotential.h"
#include "potentials/HCorePotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/DFTPotentials.h"
#include "potentials/bundles/HFPotentials.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
DensityInitialGuessCalculator::calculateInitialGuess(std::shared_ptr<SystemController> systemController) {
  /*
   * Perform the density guess
   * The guessed density is not needed anymore once the orbitals are created
   */
  auto guessDensity = calculateInitialDensity(systemController);

  auto guessOrbitals = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(
      guessDensity->getBasisController(), systemController->getNCoreElectrons() / 2);
  guessOrbitals->setCanOrthTh(systemController->getSettings().scf.canOrthThreshold);
  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> elecStruct(
      new ElectronicStructure<Options::SCF_MODES::RESTRICTED>(
          guessOrbitals, systemController->getOneElectronIntegralController(),
          systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>()));
  auto dMatController = elecStruct->getDensityMatrixController();
  double degeneracyThreshold = systemController->getSettings().scf.degeneracyThreshold;
  dMatController->setDegeneracyThreshold(degeneracyThreshold);
  auto occupations = dMatController->getOccupations(true);
  dMatController->setDensityMatrix(*guessDensity);
  dMatController->attachOrbitals(guessOrbitals, occupations, false);

  /*
   * Create Fock matrix from guessdensity
   */
  std::shared_ptr<PotentialBundle<Options::SCF_MODES::RESTRICTED>> pot;
  std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED>> pcm(
      new ZeroPotential<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));

  if (systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    std::shared_ptr<HCorePotential<Options::SCF_MODES::RESTRICTED>> hcore(
        new HCorePotential<Options::SCF_MODES::RESTRICTED>(systemController->getSharedPtr()));
    std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED>> hf(new ERIPotential<Options::SCF_MODES::RESTRICTED>(
        systemController->getSharedPtr(), dMatController, 1.0, systemController->getSettings().basis.integralThreshold,
        systemController->getSettings().basis.integralIncrementThresholdStart,
        systemController->getSettings().basis.integralIncrementThresholdEnd,
        systemController->getSettings().basis.incrementalSteps, 0.0, 0.3, false));
    pot = std::make_shared<HFPotentials<Options::SCF_MODES::RESTRICTED>>(hcore, hf, pcm, systemController->getGeometry());
  }
  else {
    // Hcore
    std::shared_ptr<HCorePotential<Options::SCF_MODES::RESTRICTED>> hcore(
        new HCorePotential<Options::SCF_MODES::RESTRICTED>(systemController->getSharedPtr()));
    // XC Func
    auto functional = systemController->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(systemController->getSettings().customFunc)
                          : resolveFunctional(systemController->getSettings().dft.functional);
    std::shared_ptr<FuncPotential<Options::SCF_MODES::RESTRICTED>> Vxc(new FuncPotential<Options::SCF_MODES::RESTRICTED>(
        systemController, dMatController, systemController->getGridController(), functional));
    // J
    std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED>> J;
    double thresh = systemController->getSettings().basis.integralThreshold;
    J = std::shared_ptr<Potential<Options::SCF_MODES::RESTRICTED>>(new ERIPotential<Options::SCF_MODES::RESTRICTED>(
        systemController->getSharedPtr(), dMatController, functional.getHfExchangeRatio(), thresh,
        systemController->getSettings().basis.integralIncrementThresholdStart,
        systemController->getSettings().basis.integralIncrementThresholdEnd,
        systemController->getSettings().basis.incrementalSteps, 0.0, 0.3, false));

    // Bundle
    pot = std::make_shared<DFTPotentials<Options::SCF_MODES::RESTRICTED>>(
        hcore, J, Vxc, pcm, systemController->getGeometry(), dMatController,
        systemController->getSettings().basis.integralThreshold);
  }
  /*
   * Create orbitals from Fock matrix
   */
  guessOrbitals->updateOrbitals(pot->getFockMatrix(*guessDensity, elecStruct->getEnergyComponentController()),
                                systemController->getOneElectronIntegralController());
  return elecStruct;
}

} /* namespace Serenity */
