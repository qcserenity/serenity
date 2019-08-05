/**
 * @file CoupledClusterTask.cpp
 *
 * @date Apr 3, 2016
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/CoupledClusterTask.h"
/* Include Serenity Internal Headers */
#include "postHF/CC/CCSD.h"
#include "postHF/CC/CCSD_T.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "system/SystemController.h"


namespace Serenity {

CoupledClusterTask::CoupledClusterTask(std::shared_ptr<SystemController> systemController) :
            _systemController(systemController){
  assert(_systemController);
}

void CoupledClusterTask::run(){

  // runs scf if need be
  _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  // Building a CCSD(T) object even if CCSD is requested,
  //  because it can do both and does not allocate any extra memory.
  // Change this if it does not make sense anymore.
  CCSD_T ccsd(_systemController, settings.normThreshold, settings.maxCycles);

  // calculate CCSD energy correction (or any other in the future)
  auto result = ccsd.calculateElectronicEnergyCorrections();

  // add energies to the EnergyController
  auto energyController = _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::MP2_CORRECTION,result.first));
  energyController->addOrReplaceComponent(std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION,result.second));
  double triples = 0.0;
  if (settings.level == Options::CC_LEVEL::CCSD_T){
    printf("\n");
    printSmallCaption("Triples Correction Calculation");
    // calculate triples energy correction
    triples = ccsd.calculateTripplesCorrection();
    // add energy to the EnergyController
    energyController->addOrReplaceComponent(std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION,triples));
    printf("%4s %10s \n","","Done.");
  }

  // get scf energy
  double energy = _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
          getEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY);

  if (energyController->checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)){
    energy = _systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
              getEnergy(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT);
  }

  /*
   * Print final results
   */
  printSectionTitle("CC Results");
  print((std::string)"Old Energy                "+energy+" Hartree\n");
  print((std::string)"MP2 Energy Correction     "+result.first+" Hartree\n");
  print((std::string)"Total MP2 Energy          "+(energy+result.first)+" Hartree\n");
  print((std::string)"CCSD Energy Correction    "+result.second+" Hartree\n");
  if (settings.level == Options::CC_LEVEL::CCSD_T){
    print((std::string)"Triples Energy Correction "+triples+" Hartree\n");
  }
  print((std::string)"Total CCSD Energy         "+(energy+result.second)+" Hartree\n");
  if (settings.level == Options::CC_LEVEL::CCSD_T){
    print((std::string)"Total CCSD(T) Energy      "+(energy+result.second+triples)+" Hartree\n");
  }
}
} /* namespace Serenity */
