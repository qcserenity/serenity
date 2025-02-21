/**
 * @file EvaluateEnergyTask.cpp
 *
 * @author Moritz Bensberg, Anja Massolle
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
#include "data/ElectronicStructure.h" //Access to the EnergyComponentController.
#include "data/OrbitalController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h" //Dispersion correction
#include "energies/EnergyComponentController.h"                      //Energy printing and adding.
#include "energies/EnergyContributions.h"                            //EnergyContributions definitions.
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h" //Kinetic integrals.
#include "io/FormattedOutput.h"                      //Captions.
#include "misc/SerenityError.h"                      //Error messages.
#include "postHF/MPn/LocalMP2.h"                     //Local MP2.
#include "postHF/MPn/MP2.h"                          //MP2.
#include "postHF/MPn/RIMP2.h"                        //RI-MP2.
#include "settings/Settings.h"
#include "settings/Settings.h"       //Settings-->HF vs DFT
#include "system/SystemController.h" //SystemController definition.
#include "tasks/FDETask.h"
#include "tasks/LocalizationTask.h" //Orbital localization for local MP2.
#include "tasks/OrthogonalizationTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemAdditionTask.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
EvaluateEnergyTask<SCFMode>::EvaluateEnergyTask(std::vector<std::shared_ptr<SystemController>> systemController,
                                                std::shared_ptr<SystemController> superSystem)
  : _systemController(systemController), _superSystem(superSystem) {
}

template<Options::SCF_MODES SCFMode>
void EvaluateEnergyTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _systemController, {_superSystem});
  if (_systemController.size() == 1) {
    printSubSectionTitle((std::string) "Energy Evaluation for System " + _systemController[0]->getSystemName());

    /*
     * KS-DFT
     */
    auto system = _systemController[0];
    auto originalXCfunc = (system->getSettings()).dft.functional;
    auto originalCustomXCfunc = system->getSettings().customFunc;
    if (settings.useDifferentXCFunc) {
      system->setXCfunctional(settings.XCfunctional);
      if (settings.customFunc.basicFunctionals.size()) {
        system->setXCfunctional(settings.customFunc);
      }
    }
    ScfTask<SCFMode> scf(_systemController[0]);
    scf.settings.skipSCF = true;
    scf.settings.mp2Type = settings.mp2Type;
    scf.settings.lcSettings = settings.lcSettings;
    scf.settings.maxResidual = settings.maxResidual;
    scf.settings.maxCycles = settings.maxCycles;
    scf.run();
    system->setXCfunctional(originalXCfunc);
    system->setXCfunctional(originalCustomXCfunc);
  }
  /*
   * sDFT evaluation
   */
  else {
    std::shared_ptr<SystemController> activeSystem;
    std::vector<std::shared_ptr<SystemController>> envSystems;
    std::vector<CompositeFunctionals::XCFUNCTIONALS> originalXCfunc;
    std::vector<CUSTOMFUNCTIONAL> originalCustomXCfunc;
    /*
     * construct Systems with the new XC func
     */

    for (unsigned int i = 0; i < _systemController.size(); i++) {
      originalXCfunc.push_back((_systemController[i]->getSettings()).dft.functional);
      originalCustomXCfunc.push_back(_systemController[i]->getSettings().customFunc);
      if (settings.useDifferentXCFunc) {
        _systemController[i]->setXCfunctional(settings.XCfunctional);
        if (settings.customFunc.basicFunctionals.size()) {
          _systemController[i]->setXCfunctional(settings.customFunc);
        }
      }

      // Update the energies of the subsystems
      ScfTask<SCFMode> scf(_systemController[i]);
      scf.settings.skipSCF = true;
      scf.run();

      if (i > 0) {
        envSystems.push_back(_systemController[i]);
      }
      else {
        activeSystem = _systemController[i];
      }
    }

    FDETask<SCFMode> fde(activeSystem, envSystems);
    fde.settings.skipSCF = true;
    fde.settings.lcSettings = this->settings.lcSettings;
    fde.settings.embedding = this->settings.embedding;
    fde.run();
    auto eCont = activeSystem->getElectronicStructure<SCFMode>()->getEnergyComponentController();

    /*
     * Evaluate the non-additive kinetic energy from orthogonalized orbitals
     */

    if (settings.evalTsOrtho) {
      double naddKin = this->calcNaddKin(_superSystem, _systemController);
      eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC, naddKin);

      /*
       * Evaluate all energy contributions from orthogonalized orbitals
       */
    }
    else if (settings.evalAllOrtho) {
      Settings settingsSuper;
      if (_superSystem == nullptr) {
        settingsSuper = _systemController[0]->getSettings();
        settingsSuper.name = _systemController[0]->getSettings().name + "Ortho";
        settingsSuper.charge = 0;
        settingsSuper.spin = 0;
        _superSystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSuper);
      }

      OrthogonalizationTask<SCFMode> orthoTask(_systemController, _superSystem);
      orthoTask.settings.orthogonalizationScheme = settings.orthogonalizationScheme;
      orthoTask.run();

      auto energyEval = EvaluateEnergyTask({_superSystem});
      energyEval.settings.XCfunctional = settings.XCfunctional;
      energyEval.run();
    }
    eCont->printAllComponents();
    for (unsigned int i = 0; i < _systemController.size(); i++) {
      _systemController[i]->setXCfunctional(originalXCfunc[i]);
      _systemController[i]->setXCfunctional(originalCustomXCfunc[i]);
    }
  }
}

template<Options::SCF_MODES SCFMode>
double EvaluateEnergyTask<SCFMode>::calcNaddKin(std::shared_ptr<SystemController> supersystem,
                                                std::vector<std::shared_ptr<SystemController>> subsystems) {
  Settings settingsSuper;
  if (supersystem == nullptr) {
    settingsSuper = _systemController[0]->getSettings();
    settingsSuper.name = _systemController[0]->getSettings().name + "Ortho";
    settingsSuper.charge = 0;
    settingsSuper.spin = 0;
    supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), settingsSuper);
  }

  if (settings.orthogonalizationScheme != Options::ORTHOGONALIZATION_ALGORITHMS::NONE) {
    OrthogonalizationTask<SCFMode> orthoTask(subsystems, supersystem);
    orthoTask.settings.orthogonalizationScheme = settings.orthogonalizationScheme;
    orthoTask.run();
  }
  else {
    SystemAdditionTask<SCFMode> additionTask(supersystem, _systemController);
    additionTask.settings.addOccupiedOrbitals = true;
    additionTask.run();
  }

  Eigen::MatrixXd sAO = supersystem->getOneElectronIntegralController()->getOverlapIntegrals();
  auto coeffMatrix = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  auto kin = supersystem->getOneElectronIntegralController()->getKinIntegrals();
  auto P = supersystem->getElectronicStructure<SCFMode>()->getDensityMatrix();

  double naddKinE = 0.0;
  for_spin(coeffMatrix, P) {
    if (settings.orthogonalizationScheme == Options::ORTHOGONALIZATION_ALGORITHMS::NONE) {
      Eigen::MatrixXd sMO = (coeffMatrix_spin).transpose() * sAO * coeffMatrix_spin;
      P_spin = (coeffMatrix_spin) * (sMO.completeOrthogonalDecomposition().pseudoInverse()).eval() *
               (coeffMatrix_spin).transpose();
    }
    naddKinE += (P_spin).cwiseProduct(kin).sum();
  };
  for (unsigned int i = 0; i < subsystems.size(); i++) {
    auto subsystemKinInt = subsystems[i]->getOneElectronIntegralController()->getKinIntegrals();
    auto subsystemDens = subsystems[i]->getElectronicStructure<SCFMode>()->getDensityMatrix();
    for_spin(subsystemDens) {
      naddKinE -= (subsystemDens_spin).cwiseProduct(subsystemKinInt).sum();
    };
  }
  return naddKinE;
} /*calcNaddKin*/

template class EvaluateEnergyTask<Options::SCF_MODES::RESTRICTED>;
template class EvaluateEnergyTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
