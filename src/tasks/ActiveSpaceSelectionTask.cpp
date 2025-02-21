/**
 * @file ActiveSpaceSelectionTask.cpp
 *
 * @date Sep 11, 2018
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
#include "tasks/ActiveSpaceSelectionTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"        //Subsystem construction.
#include "settings/Settings.h"        //Subsystem construction.
#include "system/SystemController.h"  //Subsystem construction.
#include "tasks/GeneralizedDOSTask.h" //DOS procedure.
#include "tasks/ScfTask.h"            //Initial SCF.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ActiveSpaceSelectionTask<SCFMode>::ActiveSpaceSelectionTask(std::vector<std::shared_ptr<SystemController>> supersystems,
                                                            std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                            std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _supersystems(supersystems), _activeSystems(activeSystems), _environmentSystems(environmentSystems) {
  if (_supersystems.size() <= 1)
    throw SerenityError("ERROR: Comparison between systems is only sensible if there are more than one!");
}
template<Options::SCF_MODES SCFMode>
std::shared_ptr<SystemController>
ActiveSpaceSelectionTask<SCFMode>::initializeSubsystem(std::shared_ptr<SystemController> supersystemController,
                                                       std::string namePostfix) {
  Settings subSettings = supersystemController->getSettings();
  subSettings.charge = 0;
  subSettings.spin = 0;
  subSettings.name = subSettings.name + namePostfix;
  subSettings.path =
      subSettings.path.substr(0, subSettings.path.size() - (supersystemController->getSystemName().size() + 1));
  subSettings.identifier = subSettings.identifier + namePostfix;
  return std::make_shared<SystemController>(std::make_shared<Geometry>(), subSettings);
}
template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<SystemController>> ActiveSpaceSelectionTask<SCFMode>::sortOrCreateAllSubsystemController() {
  std::vector<std::shared_ptr<SystemController>> allSubsystems;
  for (unsigned int iSuper = 0; iSuper < _supersystems.size(); ++iSuper) {
    std::shared_ptr<SystemController> activeSystemController;
    std::shared_ptr<SystemController> environmentSystemController;
    auto system = _supersystems[iSuper];
    // Create active system if needed.
    if (iSuper + 1 > _activeSystems.size()) {
      _activeSystems.push_back(initializeSubsystem(system, "_Act"));
    }
    // Create environment system if needed.
    if (iSuper + 1 > _environmentSystems.size()) {
      _environmentSystems.push_back(initializeSubsystem(system, "_Env"));
    }
    // Add systems to vector.
    allSubsystems.push_back(_activeSystems[iSuper]);
    allSubsystems.push_back(_environmentSystems[iSuper]);
    // Save the system pairs if required.
    if (keepSystemPairs) {
      _systemPairs.push_back(std::pair<std::shared_ptr<SystemController>, std::shared_ptr<SystemController>>(
          _activeSystems[iSuper], _environmentSystems[iSuper]));
    }
  } // for iSuper
  return allSubsystems;
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::run() {
  for (auto& sys : _supersystems) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The supersystem '" + sys->getSystemName() + "' needs the same SCFMode as the other systems.");
  }
  for (auto& sys : _activeSystems) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The active system '" + sys->getSystemName() + "' needs the same SCFMode as the other systems.");
  }
  for (auto& sys : _environmentSystems) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The environment system '" + sys->getSystemName() +
                          "' needs the same SCFMode as the other systems.");
  }
  // Run SCF, alignment and orbital localization.
  if (!settings.skipLocalization)
    prepareOrbitals();
  // Create or map all subsystem controller.
  auto allSubsystemController = sortOrCreateAllSubsystemController();
  // Run the generalized DOS.
  GeneralizedDOSTask<SCFMode> gdos(_supersystems, allSubsystemController);
  gdos.settings.similarityKinEnergyThreshold = {settings.similarityKinEnergyThreshold};
  gdos.settings.similarityLocThreshold = {settings.similarityLocThreshold};
  gdos.settings.prioFirst = true;
  gdos.settings.localizationThreshold = settings.localizationThreshold;
  gdos.settings.populationAlgorithm = settings.populationAlgorithm;
  gdos.settings.usePiBias = settings.usePiBias;
  gdos.settings.biasThreshold = settings.biasThreshold;
  gdos.settings.biasAverage = settings.biasAverage;
  gdos.run();
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::prepareOrbitals() {
  for (auto& sys : _supersystems) {
    if (!settings.load) {
      ScfTask<SCFMode> scfTask(sys);
      scfTask.run();
    }
    else {
      if (!sys->template hasElectronicStructure<SCFMode>())
        throw SerenityError(
            (std::string) "ERROR: No electronic structure available. However load=true was set! System " +
            sys->getSystemName());
    } // else if !settings.load
  }   // for sys
  {
    LocalizationTask locTask(_supersystems[0]);
    locTask.settings = settings.loc;
    locTask.run();
  }

  if (settings.alignPiOrbitals)
    alignPiOrbitals();

  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    LocalizationTask locTask(_supersystems[i]);
    locTask.settings = settings.loc;
    locTask.run();
  }
}

template<Options::SCF_MODES SCFMode>
void ActiveSpaceSelectionTask<SCFMode>::alignPiOrbitals() {
  for (unsigned int i = 1; i < _supersystems.size(); ++i) {
    LocalizationTask alignTask(_supersystems[i], {_supersystems[0]});
    alignTask.settings = settings.loc;
    alignTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN;
    alignTask.run();
  }
}

template class ActiveSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>;
template class ActiveSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
