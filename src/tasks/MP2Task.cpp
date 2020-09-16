/**
 * @file   MP2Task.cpp
 *
 * @date   Jul 14, 2014
 * @author Jan Unsleber
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
#include "tasks/MP2Task.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "io/FormattedOutput.h"
#include "postHF/MPn/LocalMP2.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <iomanip> //setw(...) for ostream

namespace Serenity {

template<Options::SCF_MODES SCFMode>
MP2Task<SCFMode>::MP2Task(std::shared_ptr<SystemController> systemController,
                          std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _systemController(systemController), _environmentSystems(environmentSystems) {
  assert(_systemController);
}

template<Options::SCF_MODES SCFMode>
void MP2Task<SCFMode>::run() {
  Eigen::VectorXd mp2EnergyCorrections(1);
  if (settings.mp2Type == Options::MP2_TYPES::RI) {
    RIMP2<SCFMode> rimp2(_systemController);
    mp2EnergyCorrections << rimp2.calculateCorrection();
  }
  else if (settings.mp2Type == Options::MP2_TYPES::AO) {
    if (SCFMode == UNRESTRICTED)
      throw SerenityError("MP2 is not available for unrestricted systems, please use RI-MP2.");
    MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_systemController);
    mp2EnergyCorrections << mp2EnergyCorrector.calculateElectronicEnergy();
  }
  else if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    if (SCFMode == UNRESTRICTED)
      throw SerenityError("Local-MP2 is not available for unrestricted systems, please use RI-MP2.");
    settings.lcSettings.method = PNO_METHOD::DLPNO_MP2;
    auto localCorrelationController =
        std::make_shared<LocalCorrelationController>(_systemController, settings.lcSettings, _environmentSystems);
    LocalMP2 mp2EnergyCorrector(localCorrelationController);
    mp2EnergyCorrector.settings.maxResidual = settings.maxResidual;
    mp2EnergyCorrector.settings.maxCycles = settings.maxCycles;
    mp2EnergyCorrections = mp2EnergyCorrector.calculateEnergyCorrection();
  }

  // add energy to the EnergyController
  auto energyController = _systemController->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::MP2_CORRECTION, mp2EnergyCorrections.sum()));

  /*
   * Print final results
   */
  if (settings.mp2Type == Options::MP2_TYPES::RI) {
    printSectionTitle("(RI-)MP2 Results");
  }
  else if (settings.mp2Type == Options::MP2_TYPES::AO) {
    printSectionTitle("MP2 Results");
  }
  else if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    printSectionTitle("(Local-)MP2 Results");
  }
  double energy = _systemController->getElectronicStructure<SCFMode>()->getEnergy() - mp2EnergyCorrections.sum();
  std::cout << std::fixed;
  std::cout << "Old Energy                        " << setw(20) << energy << " Hartree   ";
  std::cout << setw(20) << energy * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
  if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    std::cout << "Local MP2 pair energies           " << setw(20) << mp2EnergyCorrections(0) << " Hartree   ";
    std::cout << setw(20) << mp2EnergyCorrections(0) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "Dipole app. contribution          " << setw(20) << mp2EnergyCorrections(1) << " Hartree   ";
    std::cout << setw(20) << mp2EnergyCorrections(1) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "PNO truncation correction         " << setw(20) << mp2EnergyCorrections(2) << " Hartree   ";
    std::cout << setw(20) << mp2EnergyCorrections(2) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "Unscaled total Local-MP2 energy   " << setw(20) << mp2EnergyCorrections.sum() << " Hartree   ";
    std::cout << setw(20) << mp2EnergyCorrections.sum() * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
  }
  else {
    std::cout << "Canonical MP2 energy correction   " << setw(20) << mp2EnergyCorrections(0) << " Hartree   ";
    std::cout << setw(20) << mp2EnergyCorrections(0) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
  }
  double totalEnergy = _systemController->getElectronicStructure<SCFMode>()->getEnergy();
  std::cout << "Total Energy                      " << setw(20) << totalEnergy << " Hartree   ";
  std::cout << setw(20) << totalEnergy * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
  std::cout << std::scientific;
  return;
}

template class MP2Task<Options::SCF_MODES::RESTRICTED>;
template class MP2Task<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
