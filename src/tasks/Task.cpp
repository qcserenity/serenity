/**
 * @file   Task.cpp
 *
 * @date   Feb. 10, 2025
 * @author Niklas Goellmann
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
#include "tasks/Task.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"

namespace Serenity {

void Task::avoidMixedSCFModes(Options::SCF_MODES runMode, std::vector<std::shared_ptr<SystemController>> activeSystems,
                              std::vector<std::shared_ptr<SystemController>> superSystems) {
  if (runMode == Options::SCF_MODES::UNRESTRICTED) {
    for (auto system : activeSystems) {
      if (!system)
        continue;
      if (system->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        system->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
        if (system->hasElectronicStructure<RESTRICTED>())
          system->setElectronicStructure(
              std::make_shared<ElectronicStructure<UNRESTRICTED>>(system->getElectronicStructure<RESTRICTED>()));
        WarningTracker::printWarning("Changing SCFMode of system '" + system->getSystemName() + "' to unrestricted.", true);
      }
    }
    for (auto system : superSystems) {
      if (!system)
        continue;
      if (system->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        system->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
        if (system->hasElectronicStructure<RESTRICTED>())
          system->setElectronicStructure(
              std::make_shared<ElectronicStructure<UNRESTRICTED>>(system->getElectronicStructure<RESTRICTED>()));
        WarningTracker::printWarning("Changing SCFMode of system '" + system->getSystemName() + "' to unrestricted.", true);
      }
    }
  }
}

void Task::avoidMixedSCFModes(Options::SCF_MODES runMode, std::vector<std::shared_ptr<SystemController>> activeSystems) {
  avoidMixedSCFModes(runMode, activeSystems, {});
}

void Task::avoidMixedSCFModes(Options::SCF_MODES runMode, std::shared_ptr<SystemController> activeSystem) {
  avoidMixedSCFModes(runMode, {activeSystem}, {});
}

} /* namespace Serenity */
