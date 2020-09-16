/**
 * @file CoupledClusterTask.cpp
 *
 * @date Apr 3, 2016
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
#include "tasks/CoupledClusterTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"           //Get and add energies.
#include "energies/EnergyComponentController.h" //Add energies to the system.
#include "energies/EnergyContributions.h"       //Get/add specific energy contributions.
#include "io/FormattedOutput.h"                 //Captions.
#include "misc/WarningTracker.h"                //Warnings for incorrect energy handling.
#include "postHF/CC/CCSD.h"                     //Canonical coupled cluster.
#include "postHF/CC/CCSD_T.h"                   //Canonical coupled cluster.
#include "postHF/CC/DLPNO_CCSD.h"               //Local coupled cluster.
#include "postHF/CC/DLPNO_CCSD_T0.h"            //Local coupled cluster
#include "system/SystemController.h"            //Access to the electronic structure of the system.
/* Include Std and External Headers */
#include <iomanip> //setw(...) for ostream

namespace Serenity {

CoupledClusterTask::CoupledClusterTask(std::shared_ptr<SystemController> systemController,
                                       std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _systemController(systemController), _environmentSystems(environmentSystems) {
  assert(_systemController);
}

Eigen::VectorXd CoupledClusterTask::localCalculation() {
  if (settings.level == Options::CC_LEVEL::DLPNO_CCSD_T0)
    settings.lcSettings.method = PNO_METHOD::DLPNO_CCSD_T0;
  if (settings.level == Options::CC_LEVEL::DLPNO_CCSD)
    settings.lcSettings.method = PNO_METHOD::DLPNO_CCSD;
  auto localCorrelationController =
      std::make_shared<LocalCorrelationController>(_systemController, settings.lcSettings, _environmentSystems);
  DLPNO_CCSD dlpnoCCSD(localCorrelationController, settings.normThreshold, settings.maxCycles,
                       settings.level == Options::CC_LEVEL::DLPNO_CCSD_T0);
  Eigen::VectorXd energies = dlpnoCCSD.calculateElectronicEnergyCorrections();
  auto energyController =
      _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION, energies.sum()));
  if (settings.level == Options::CC_LEVEL::DLPNO_CCSD_T0) {
    energies.conservativeResize(energies.size() + 1);
    energies[5] = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
    energyController->addOrReplaceComponent(
        std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION, energies[5]));
  }
  return energies;
}

Eigen::VectorXd CoupledClusterTask::canonicalCalculation() {
  // Building a CCSD(T) object even if CCSD is requested,
  // because it can do both and does not allocate any extra memory.
  // Change this if it does not make sense anymore.
  CCSD_T ccsd(_systemController, settings.normThreshold, settings.maxCycles);

  // calculate CCSD energy correction (or any other in the future)
  auto result = ccsd.calculateElectronicEnergyCorrections();

  // add energies to the EnergyController
  auto energyController =
      _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::MP2_CORRECTION, result.first));
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION, result.second));
  double triples = 0.0;
  if (settings.level == Options::CC_LEVEL::CCSD_T) {
    printf("\n");
    printSmallCaption("Triples Correction Calculation");
    // calculate triples energy correction
    triples = ccsd.calculateTripplesCorrection();
    // add energy to the EnergyController
    energyController->addOrReplaceComponent(
        std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION, triples));
    printf("%4s %10s \n", "", "Done.");
  }
  Eigen::VectorXd energies(3);
  energies << result.first, result.second, triples;
  return energies;
}

void CoupledClusterTask::run() {
  bool canonical = settings.level == Options::CC_LEVEL::CCSD || settings.level == Options::CC_LEVEL::CCSD_T;
  bool dlpno = !canonical;
  // runs scf if need be
  _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  Eigen::VectorXd energies;
  if (canonical)
    energies = canonicalCalculation();
  if (dlpno)
    energies = localCalculation();

  // get scf energy
  auto es = _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  double energy = es->getEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY);
  if (_environmentSystems.size() > 0 || settings.level == Options::CC_LEVEL::CCSD ||
      settings.level == Options::CC_LEVEL::CCSD_T) {
    if (es->checkEnergy(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)) {
      energy = es->getEnergy(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT);
    }
    else if (es->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)) {
      energy = es->getEnergy(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
    }
    else {
      if (_environmentSystems.size() > 0)
        WarningTracker::printWarning(
            (std::string) "Warning: This calculation was set up as an embedded coupled cluster calculation.\n" +
                "         However, the active system has not the full information about all energy\n" +
                "         contributions that originate from the embedding. Handle the following energies\n" +
                "         with care. Typical errors are:\n" +
                "            * Did you specify an environment system that is not needed?\n" +
                "            * Did you forget to calculate the energy of the frozen environment by not setting" +
                "                calculateEnvironmentEnergy = true in a previous FDETask?" +
                "            * Is the method used for the active system set to DFT?",
            true);
    }
  }

  /*
   * Print final results
   */
  std::cout << std::fixed;
  printSectionTitle("CC Results");
  std::cout << "Old Energy                        " << setw(20) << energy << " Hartree   ";
  std::cout << setw(20) << energy * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  if (canonical) {
    std::cout << "MP2 Energy Correction             " << setw(20) << energies(0) << " Hartree   ";
    std::cout << setw(20) << energies(0) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "Total MP2 Energy                  " << setw(20) << (energy + energies(2)) << " Hartree   ";
    std::cout << setw(20) << (energy + energies(2)) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "CCSD Energy Correction            " << setw(20) << energies(1) << " Hartree   ";
    std::cout << setw(20) << energies(1) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    if (settings.level == Options::CC_LEVEL::CCSD_T) {
      std::cout << "Triples Energy Correction         " << setw(20) << energies(2) << " Hartree   ";
      std::cout << setw(20) << energies(2) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    if (settings.level == Options::CC_LEVEL::CCSD_T) {
      std::cout << "Total CCSD(T) Energy              " << setw(20) << energies(1) + energies(2) + energy << " Hartree   ";
      std::cout << setw(20) << (energies(1) + energies(2) + energy) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    }
    else {
      std::cout << "Total CCSD Energy                 " << setw(20) << energies(1) + energy << " Hartree   ";
      std::cout << setw(20) << (energies(1) + energy) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    }
  }
  else {
    std::cout << "Local CCSD pair energies          " << setw(20) << energies(0) << " Hartree   ";
    std::cout << setw(20) << energies(0) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "SC-MP2 pair energies              " << setw(20) << energies(2) << " Hartree   ";
    std::cout << setw(20) << energies(2) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "Dipole app. contribution          " << setw(20) << energies(3) << " Hartree   ";
    std::cout << setw(20) << energies(3) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    std::cout << "PNO truncation correction         " << setw(20) << energies(4) << " Hartree   ";
    std::cout << setw(20) << energies(4) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    if (settings.level == Options::CC_LEVEL::DLPNO_CCSD) {
      std::cout << "Local-CCSD Energy Correction      " << setw(20) << energies.segment(0, energies.size()).sum()
                << " Hartree   ";
      std::cout << setw(20) << energies.segment(0, 4).sum() * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
      std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "Total Local-CCSD Energy           " << setw(20) << energies.sum() + energy << " Hartree   ";
      std::cout << setw(20) << (energies.sum() + energy) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    }
    else {
      std::cout << "Local-CCSD Energy Correction      " << setw(20) << energies.segment(0, energies.size() - 1).sum()
                << " Hartree   ";
      std::cout << setw(20) << energies.segment(0, 4).sum() * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
      std::cout << "Semi-Canonical Triples Correction " << setw(20) << energies(5) << " Hartree   ";
      std::cout << setw(20) << energies(5) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
      std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "Total Local-CCSD(T0) Energy       " << setw(20) << energies.sum() + energy << " Hartree   ";
      std::cout << setw(20) << (energies.sum() + energy) * HARTREE_TO_KJ_PER_MOL << " kJ/mol" << std::endl;
    }
  }
  std::cout << std::scientific << std::endl;
}
} /* namespace Serenity */
