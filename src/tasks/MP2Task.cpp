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
#include "energies/EnergyComponentController.h" //Assign and get energy contributions.
#include "energies/EnergyContributions.h"       //Assign and get energy contributions.
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h" //Warnings for incorrect energy handling.
#include "postHF/MPn/LTMP2.h"    //Laplace-Transform SOS RI MP2.
#include "postHF/MPn/LocalMP2.h" //Local MP2.
#include "postHF/MPn/MP2.h"      //Canonical MP2.
#include "postHF/MPn/RIMP2.h"    //RI MP2.
#include "settings/Settings.h"
#include "system/SystemController.h" //Access to electronic structure.
/* Include Std and External Headers */
#include <iomanip> //std::setw(...) for ostream

namespace Serenity {

template<Options::SCF_MODES SCFMode>
MP2Task<SCFMode>::MP2Task(std::shared_ptr<SystemController> systemController,
                          std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _systemController(systemController), _environmentSystems(environmentSystems) {
  assert(_systemController);
}

template<Options::SCF_MODES SCFMode>
void MP2Task<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _systemController);
  Eigen::VectorXd mp2EnergyCorrections(1);
  if (settings.sss != 1.0 || settings.oss != 1.0) {
    OutputControl::nOut << "Custom spin-component scaling:" << std::endl;
    OutputControl::nOut << " Same-spin     : " << settings.sss << std::endl;
    OutputControl::nOut << " Opposite-spin : " << settings.oss << std::endl;
  }
  if (settings.ltconv != 0) {
    printBigCaption("Laplace-Transform algorithms");
    printf(" Threshold     : %-6.1e\n", settings.ltconv);
    printf(" Opposite-spin : %-6.3f\n\n", settings.oss);
    if (settings.mp2Type != Options::MP2_TYPES::LT) {
      WarningTracker::printWarning("You have set the ltconv threshold but did not request LT as the MP2 type.", 1);
    };
  }
  if (settings.mp2Type == Options::MP2_TYPES::DF) {
    RIMP2<SCFMode> rimp2(_systemController, settings.sss, settings.oss);
    mp2EnergyCorrections << rimp2.calculateCorrection();
    if (settings.unrelaxedDensity) {
      auto densityCorrection = rimp2.calculateDensityCorrection();
      auto densityController = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
      auto density = densityController->getDensityMatrix();
      densityController->setDensityMatrix(density + densityCorrection);
    }
  }
  else if (settings.mp2Type == Options::MP2_TYPES::LT) {
    LTMP2<SCFMode> ltmp2(_systemController, settings.oss, settings.ltconv);
    mp2EnergyCorrections << ltmp2.calculateCorrection();
  }
  else if (settings.mp2Type == Options::MP2_TYPES::AO) {
    if (settings.sss != 1.0 || settings.sss != 1.0) {
      WarningTracker::printWarning("Custom spin scaling is ignored.", 1);
    }
    if (SCFMode == UNRESTRICTED) {
      throw SerenityError("MP2 is not available for unrestricted systems, please use RI-MP2.");
    }
    MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_systemController);
    mp2EnergyCorrections << mp2EnergyCorrector.calculateElectronicEnergy();
    if (settings.unrelaxedDensity) {
      throw SerenityError("Unrelaxed density is not available for canonical MP2, please use RI-MP2.");
    }
  }
  else if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    if (SCFMode == UNRESTRICTED) {
      throw SerenityError("Local-MP2 is not available for unrestricted systems, please use RI-MP2.");
    }
    // If SC-MP2 is not requested explicitly, set the DLPNO method to DLPNO-MP2.
    if (settings.lcSettings.method != Options::PNO_METHOD::SC_MP2) {
      settings.lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
    }
    auto localCorrelationController =
        std::make_shared<LocalCorrelationController>(_systemController, settings.lcSettings, _environmentSystems);
    LocalMP2 mp2EnergyCorrector(localCorrelationController);
    mp2EnergyCorrector.settings.ssScaling = settings.sss;
    mp2EnergyCorrector.settings.osScaling = settings.oss;
    mp2EnergyCorrector.settings.maxResidual = settings.maxResidual;
    mp2EnergyCorrector.settings.maxCycles = settings.maxCycles;
    mp2EnergyCorrections = mp2EnergyCorrector.calculateEnergyCorrection();
    if (settings.unrelaxedDensity) {
      auto densityCorrection = mp2EnergyCorrector.calculateDensityCorrection();
      auto densityController = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController();
      auto density = densityController->getDensityMatrix();
      densityController->setDensityMatrix(density + densityCorrection);
    }
    if (this->settings.writePairEnergies)
      localCorrelationController->writePairEnergies("MP2");
  }

  // add energy to the EnergyController
  auto energyController = _systemController->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::MP2_CORRECTION, mp2EnergyCorrections.sum()));

  /*
   * Print final results
   */
  if (settings.mp2Type == Options::MP2_TYPES::DF) {
    printSectionTitle("(RI-)MP2 Results");
  }
  if (settings.mp2Type == Options::MP2_TYPES::LT) {
    printSectionTitle("(LT-SOS-)MP2 Results");
  }
  else if (settings.mp2Type == Options::MP2_TYPES::AO) {
    printSectionTitle("MP2 Results");
  }
  else if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    printSectionTitle("(Local-)MP2 Results");
  }
  // get scf energy
  auto es = _systemController->getElectronicStructure<SCFMode>();
  double energy = 0.0;
  if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    energy = es->getEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY);
  }
  else if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    WarningTracker::printWarning("Warning: Orbitals for the MP2 Task were calculated using DFT - A Görling-Levy "
                                 "Kohn-Sham perturbation theory calculation was performed (based on the KS orbitals). "
                                 "If you want conventional MP2 select \"method HF\" in the system settings",
                                 true);
    energy = es->getEnergy(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
  }
  else
    throw SerenityError("No fitting electronic structure method selected for MP2 Task");

  if (_environmentSystems.size() > 0 || settings.mp2Type != Options::MP2_TYPES::LOCAL) {
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

  unsigned int lineLength = 80;
  unsigned int indent = 2;
  unsigned int numberLength = lineLength - 8 - 34 - 2 * indent;
  std::cout << std::fixed;
  std::cout << std::string(indent, ' ') << "Old Energy                        " << std::setw(numberLength) << energy
            << " Hartree" << std::endl;
  if (settings.mp2Type == Options::MP2_TYPES::LOCAL) {
    std::cout << std::string(indent, ' ') << "Local MP2 pair energies           " << std::setw(numberLength)
              << mp2EnergyCorrections(0) << " Hartree" << std::endl;
    std::cout << std::string(indent, ' ') << "Dipole app. contribution          " << std::setw(numberLength)
              << mp2EnergyCorrections(1) << " Hartree" << std::endl;
    std::cout << std::string(indent, ' ') << "PNO truncation correction         " << std::setw(numberLength)
              << mp2EnergyCorrections(2) << " Hartree" << std::endl;
    std::cout << std::string(indent, ' ') << "Local-MP2 energy                  " << std::setw(numberLength)
              << mp2EnergyCorrections.sum() << " Hartree" << std::endl;
  }
  else {
    std::cout << std::string(indent, ' ') << "Canonical MP2 energy correction   " << std::setw(numberLength)
              << mp2EnergyCorrections(0) << " Hartree" << std::endl;
  }
  std::cout << std::string(lineLength, '-') << std::endl;
  std::cout << std::string(indent, ' ') << "Total Energy                      " << std::setw(numberLength)
            << energy + mp2EnergyCorrections.sum() << " Hartree" << std::endl;
  std::cout << std::scientific;
  return;
}

template class MP2Task<Options::SCF_MODES::RESTRICTED>;
template class MP2Task<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
