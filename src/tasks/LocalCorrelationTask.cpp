/**
 * @file LocalCorrelationTask.cpp
 *
 * @date Apr 4, 2021
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
#include "tasks/LocalCorrelationTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "system/SystemController.h"
#include "tasks/CoupledClusterTask.h"
#include "tasks/MP2Task.h"

namespace Serenity {
LocalCorrelationTask::LocalCorrelationTask(std::shared_ptr<SystemController> system,
                                           std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _system(system), _environmentSystems(environmentSystems) {
}

LocalCorrelationTask::~LocalCorrelationTask() = default;

void LocalCorrelationTask::run() {
  printSectionTitle("Local Correlation Calculation");
  std::string pnoMethod, locMethod, pnoSettings;
  Options::resolve<Options::PNO_METHOD>(pnoMethod, this->settings.lcSettings.method);
  Options::resolve<Options::ORBITAL_LOCALIZATION_ALGORITHMS>(locMethod, this->settings.loc.locType);
  Options::resolve<Options::PNO_SETTINGS>(pnoSettings, this->settings.lcSettings.pnoSettings);
  OutputControl::nOut << "  Correlation method:  " << pnoMethod << std::endl;
  OutputControl::nOut << "  Localization method: " << locMethod << std::endl;
  OutputControl::nOut << "  DLPNO-Thresholds:    " << pnoSettings << std::endl;

  LocalizationTask locTask(_system);
  locTask.settings = settings.loc;
  locTask.run();

  switch (settings.lcSettings.method) {
    case Options::PNO_METHOD::NONE:
      OutputControl::mOut << "No correlation method selected exiting." << std::endl;
      break;
    case Options::PNO_METHOD::SC_MP2:
    case Options::PNO_METHOD::DLPNO_MP2: {
      MP2Task<RESTRICTED> mp2Task(_system, _environmentSystems);
      mp2Task.settings.lcSettings = this->settings.lcSettings;
      mp2Task.settings.mp2Type = Options::MP2_TYPES::LOCAL;
      mp2Task.settings.maxCycles = this->settings.maxCycles;
      mp2Task.settings.maxResidual = this->settings.normThreshold;
      mp2Task.run();
      break;
    }
    case Options::PNO_METHOD::DLPNO_CCSD:
    case Options::PNO_METHOD::DLPNO_CCSD_T0:
      CoupledClusterTask ccTask(_system, _environmentSystems);
      ccTask.settings.lcSettings = settings.lcSettings;
      ccTask.settings.level = (settings.lcSettings.method == Options::PNO_METHOD::DLPNO_CCSD)
                                  ? Options::CC_LEVEL::DLPNO_CCSD
                                  : Options::CC_LEVEL::DLPNO_CCSD_T0;
      ccTask.settings.writePairEnergies = this->settings.writePairEnergies;
      ccTask.settings.maxCycles = this->settings.maxCycles;
      ccTask.settings.normThreshold = this->settings.normThreshold;
      ccTask.run();
  }
}

double LocalCorrelationTask::getCorrelationEnergy(std::shared_ptr<SystemController> activeSystemController,
                                                  Options::PNO_METHOD pnoMethod) {
  double correlationEnergy = 0.0;
  auto electronicStructure = activeSystemController->getElectronicStructure<RESTRICTED>();
  switch (pnoMethod) {
    case Options::PNO_METHOD::SC_MP2:
    case Options::PNO_METHOD::DLPNO_MP2: {
      correlationEnergy = electronicStructure->getEnergy(ENERGY_CONTRIBUTIONS::MP2_CORRECTION);
      break;
    }
    case Options::PNO_METHOD::DLPNO_CCSD: {
      correlationEnergy = electronicStructure->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
      break;
    }
    case Options::PNO_METHOD::DLPNO_CCSD_T0: {
      correlationEnergy = electronicStructure->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
      correlationEnergy += electronicStructure->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
      break;
    }
    case Options::PNO_METHOD::NONE: {
      correlationEnergy = 0.0;
      break;
    }
  }
  return correlationEnergy;
}

} /* namespace Serenity */
