/**
 * @file   PCMInteractionEnergyTask.cpp
 *
 * @date   Nov 12, 2020
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
/* Include Class Header*/
#include "tasks/PCMInteractionEnergyTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "geometry/MolecularSurfaceController.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"
#include "potentials/PCMPotential.h"
#include "solvation/Solvents.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/Eigen>

namespace Serenity {
PCMInteractionEnergyTask::PCMInteractionEnergyTask(std::shared_ptr<SystemController> systems)
  : _systemController(systems) {
}
void PCMInteractionEnergyTask::run() {
  const auto lastSCFMode = _systemController->getLastSCFMode();
  if (lastSCFMode == Options::SCF_MODES::RESTRICTED) {
    this->runSpinPolarized<RESTRICTED>();
  }
  else if (lastSCFMode == Options::SCF_MODES::UNRESTRICTED) {
    this->runSpinPolarized<UNRESTRICTED>();
  }
  else {
    throw SerenityError("ERROR: Unknown Options for SCFMode in task PCMInteractionEnergyTask.");
  }
}

template<Options::SCF_MODES SCFMode>
void PCMInteractionEnergyTask::runSpinPolarized() {
  auto pcm = std::make_shared<PCMPotential<SCFMode>>(
      this->settings.pcm, _systemController->getBasisController(), _systemController->getGeometry(),
      _systemController->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE),
      _systemController->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE_VDW),
      _systemController->getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::ACTIVE));
  auto dmat = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrix();
  auto fmat = pcm->getMatrix();
  OutputControl::mOut << "---------------------------------------------------------" << std::endl;
  Solvents::printSolventInfo(this->settings.pcm);
  const double energy = pcm->getEnergy(dmat);
  OutputControl::mOut << "---------------------------------------------------------" << std::endl;
  OutputControl::mOut << "  PCM Correction afterwards " << energy << std::endl;
  OutputControl::mOut << "---------------------------------------------------------" << std::endl;
  auto eComp = _systemController->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  eComp->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::SOLVATION_ENERGY, energy);
}

} /* namespace Serenity */
