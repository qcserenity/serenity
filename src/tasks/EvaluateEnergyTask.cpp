/**
 * @file EvaluateEnergyTask.cpp
 *
 * @author Moritz Bensberg
 * @date Mar 10, 2020
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
#include "tasks/EvaluateEnergyTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"                                //Access to the EnergyComponentController.
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h" //Dispersion correction
#include "energies/EnergyComponentController.h"                      //Energy printing and adding.
#include "energies/EnergyContributions.h"                            //EnergyContributions definitions.
#include "io/FormattedOutput.h"                                      //Captions.
#include "misc/SerenityError.h"                                      //Error messages.
#include "postHF/MPn/LocalMP2.h"                                     //Local MP2.
#include "postHF/MPn/MP2.h"                                          //MP2.
#include "postHF/MPn/RIMP2.h"                                        //RI-MP2.
#include "potentials/bundles/PotentialBundle.h"                      //Fock matrix construction for energy calculation.
#include "settings/Settings.h"                                       //Settings-->HF vs DFT
#include "system/SystemController.h"                                 //SystemController definition.
#include "tasks/LocalizationTask.h"                                  //Orbital localization for local MP2.

namespace Serenity {
template<Options::SCF_MODES SCFMode>
EvaluateEnergyTask<SCFMode>::EvaluateEnergyTask(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController) {
}

template<Options::SCF_MODES SCFMode>
void EvaluateEnergyTask<SCFMode>::run() {
  printSubSectionTitle((std::string) "Energy Evaluation for System " + _systemController->getSystemName());
  auto es = _systemController->getElectronicStructure<SCFMode>();
  es->getDensityMatrixController()->updateDensityMatrix();
  const auto p = es->getDensityMatrixController()->getDensityMatrix();

  std::shared_ptr<PotentialBundle<SCFMode>> potentials;
  const Settings& settings = _systemController->getSettings();
  if (settings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    potentials = _systemController->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  }
  else if (settings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    potentials = _systemController->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  }
  else {
    throw SerenityError("Unknown electronic structure theory!");
  }
  auto energyComponentController = es->getEnergyComponentController();
  auto f = potentials->getFockMatrix(p, energyComponentController);
  double dispersionCorrection = 0;
  if (settings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    dispersionCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
        settings.dft.dispersion, _systemController->getGeometry(), settings.dft.functional);
  }
  energyComponentController->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, dispersionCorrection);

  auto functional = resolveFunctional(settings.dft.functional);
  double MP2Correlation = 0.0;
  if (functional.isDoubleHybrid()) {
    // perform MP2 for double hybrids
    switch (this->settings.mp2Type) {
      case Options::MP2_TYPES::LOCAL: {
        LocalizationTask locTask(_systemController);
        locTask.settings.splitValenceAndCore = true;
        locTask.run();
        if (SCFMode != RESTRICTED)
          throw SerenityError("MP2 is not available for unrestricted systems, please use RI-MP2. Please set "
                              "DensityFitting to RI in the system block.");
        auto localCorrelationController =
            std::make_shared<LocalCorrelationController>(_systemController, this->settings.lcSettings);
        LocalMP2 localMP2(localCorrelationController);
        localMP2.settings.ssScaling = functional.getssScaling();
        localMP2.settings.osScaling = functional.getosScaling();
        localMP2.settings.maxCycles = this->settings.maxCycles;
        localMP2.settings.maxResidual = this->settings.maxResidual;
        MP2Correlation = localMP2.calculateEnergyCorrection().sum();
        break;
      }
      case Options::MP2_TYPES::RI: {
        RIMP2<SCFMode> rimp2(_systemController, functional.getssScaling(), functional.getosScaling());
        MP2Correlation = rimp2.calculateCorrection();
        break;
      }
      case Options::MP2_TYPES::AO: {
        if (SCFMode != RESTRICTED)
          throw SerenityError("MP2 is not available for unrestricted systems, please use RI-MP2. Please set "
                              "DensityFitting to RI in the system block.");
        MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_systemController, functional.getssScaling(),
                                                          functional.getosScaling());
        MP2Correlation = mp2EnergyCorrector.calculateElectronicEnergy();
        break;
      }
    }
    MP2Correlation *= functional.getHfCorrelRatio();
  }

  // add energy to the EnergyController
  energyComponentController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, MP2Correlation));

  energyComponentController->printAllComponents();
}

template class EvaluateEnergyTask<Options::SCF_MODES::RESTRICTED>;
template class EvaluateEnergyTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
