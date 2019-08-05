/**
 * @file   MP2Task.cpp
 *
 * @date   Jul 14, 2014
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
#include "tasks/MP2Task.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "energies/EnergyContributions.h"
#include "io/FormattedOutput.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "system/SystemController.h"


namespace Serenity {

template<Options::SCF_MODES SCFMode>
MP2Task<SCFMode>::MP2Task(std::shared_ptr<SystemController> systemController) :
    _systemController(systemController) {
  assert(_systemController);
}

template<Options::SCF_MODES SCFMode>
void MP2Task<SCFMode>::run() {
  double MP2EnergyCorrection = 0.0;
  if (settings.ri){
    RIMP2<SCFMode> rimp2(_systemController);
    MP2EnergyCorrection = rimp2.calculateCorrection();
  } else {
    assert(SCFMode==RESTRICTED && "MP2 is not available for unrestricted systems please use RI-MP2.");
    MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector( _systemController);
    MP2EnergyCorrection = mp2EnergyCorrector.calculateElectronicEnergy();
  }

  // add energy to the EnergyController
  auto energyController = _systemController->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  energyController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::MP2_CORRECTION,MP2EnergyCorrection));

  /*
   * Print final results
   */
  if (settings.ri){
    printSectionTitle("(RI-)MP2 Results");
  }else{
    printSectionTitle("MP2 Results");
  }

  printSmallCaption("Old Energy");
  double energy = _systemController->getElectronicStructure<SCFMode>()->
          getEnergy()-MP2EnergyCorrection;
  print((std::string)""+energy+" Hartree");
  print("");

  printSmallCaption("MP2 Energy Correction");
  print((std::string)""+MP2EnergyCorrection+" Hartree");
  print("");


  printSmallCaption("Corrected Energy");
  print((std::string)""+_systemController->getElectronicStructure<SCFMode>()->
          getEnergy()+" Hartree");
  print("");

  return;
}

template class MP2Task<Options::SCF_MODES::RESTRICTED>;
template class MP2Task<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
