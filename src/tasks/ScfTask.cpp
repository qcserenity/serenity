/**
 * @file   ScfTask.cpp
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#include "tasks/ScfTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "postHF/MPn/LTMP2.h"
#include "postHF/MPn/LocalMP2.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "potentials/FuncPotential.h"
#include "potentials/bundles/DFTPotentials.h"
#include "potentials/bundles/PotentialBundle.h"
#include "scf/DeltaScf.h"
#include "scf/ROHF.h"
#include "scf/SCFAnalysis.h"
#include "scf/Scf.h"
#include "settings/Settings.h"
#include "solvation/Solvents.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ScfTask<SCFMode>::ScfTask(const std::shared_ptr<SystemController> systemController)
  : _systemController(systemController) {
}
template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::calculateMP2Contribution() {
  auto& systemSettings = _systemController->getSettings();
  auto es = _systemController->getElectronicStructure<SCFMode>();
  auto energyComponentController = es->getEnergyComponentController();
  // check for double hybrid functional
  auto functional = systemSettings.customFunc.basicFunctionals.size() ? Functional(systemSettings.customFunc)
                                                                      : resolveFunctional(systemSettings.dft.functional);
  double MP2Correlation = 0.0;
  if (functional.isDoubleHybrid() && settings.calculateMP2Energy) {
    Timings::takeTime("Active System -    DH-MP2 Corr.");
    // perform MP2 for double hybrids
    switch (this->settings.mp2Type) {
      case Options::MP2_TYPES::LOCAL: {
        LocalizationTask locTask(_systemController);
        locTask.settings.splitValenceAndCore = true;
        locTask.run();
        if (SCFMode != RESTRICTED) {
          throw SerenityError("Local MP2 is not available for unrestricted systems, please use RI-MP2. Please set "
                              "<mp2Type> to <RI> in the task settings.");
        }
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
      case Options::MP2_TYPES::DF: {
        if (_systemController->getSettings().basis.densFitCorr != Options::DENS_FITS::NONE) {
          RIMP2<SCFMode> rimp2(_systemController, functional.getssScaling(), functional.getosScaling());
          MP2Correlation = rimp2.calculateCorrection();
        }
        break;
      }
      case Options::MP2_TYPES::LT: {
        LTMP2<SCFMode> ltmp2(_systemController);
        MP2Correlation = ltmp2.calculateCorrection();
        break;
      }
      case Options::MP2_TYPES::AO: {
        if (SCFMode != RESTRICTED)
          throw SerenityError("MP2 is not available for unrestricted systems, please use RI-MP2. Please set "
                              "<mp2Type> to <RI> in the task settings.");
        MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_systemController, functional.getssScaling(),
                                                          functional.getosScaling());
        MP2Correlation = mp2EnergyCorrector.calculateElectronicEnergy();
        break;
      }
    }
    MP2Correlation *= functional.getHfCorrelRatio();
    Timings::timeTaken("Active System -    DH-MP2 Corr.");
  }

  // add energy to the EnergyController
  energyComponentController->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, MP2Correlation));
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::calculateDispersionCorrection() {
  auto& systemSettings = _systemController->getSettings();
  auto es = _systemController->getElectronicStructure<SCFMode>();
  auto energyComponentController = es->getEnergyComponentController();
  auto functional = systemSettings.customFunc.basicFunctionals.size() ? Functional(systemSettings.customFunc)
                                                                      : resolveFunctional(systemSettings.dft.functional);
  auto dispCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      systemSettings.dft.dispersion, _systemController->getGeometry(), systemSettings.dft.functional);
  es->getEnergyComponentController()->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, dispCorrection);
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::finalDFTEnergyEvaluation(std::shared_ptr<SPMatrix<SCFMode>> momMatrix) {
  const Settings& systemSettings = _systemController->getSettings();
  if (systemSettings.grid.accuracy != systemSettings.grid.smallGridAccuracy) {
    auto es = _systemController->getElectronicStructure<SCFMode>();
    auto energyComponentController = es->getEnergyComponentController();
    // XC Func with default grid.
    auto functional = systemSettings.customFunc.basicFunctionals.size() ? Functional(systemSettings.customFunc)
                                                                        : resolveFunctional(systemSettings.dft.functional);
    std::shared_ptr<FuncPotential<SCFMode>> Vxc = std::make_shared<FuncPotential<SCFMode>>(
        _systemController, es->getDensityMatrixController(), _systemController->getGridController(), functional);
    // Replace functional Fock matrix contribution.
    // We know that this will be an object of type DFTPotentials. Otherwise this cast would
    // be unsafe.
    auto dftPotentials = es->getPotentialBundle();
    static_cast<DFTPotentials<SCFMode>&>(*dftPotentials).replaceFunctionalPotential(Vxc);
    // Update orbitals/Fock matrix.
    auto orbitalController = es->getMolecularOrbitals();
    auto f = dftPotentials->getFockMatrix(es->getDensityMatrix(), energyComponentController);
    if (systemSettings.scf.rohf != Options::ROHF_TYPES::NONE) {
      ROHF<SCFMode>::addConstraint(f, es, systemSettings.scf.rohf, systemSettings.scf.suhfLambda);
    }
    orbitalController->updateOrbitals(f, es->getOneElectronIntegralController(), momMatrix);
    // Set new Fock matrix of the final energy evaluation
    es->setFockMatrix(f);
    // Write everything to HDF5
    es->toHDF5(systemSettings.path + systemSettings.name, systemSettings.identifier);
  }
  // Double hybrid MP2 contribution.
  calculateMP2Contribution();
  // Dispersion correction.
  calculateDispersionCorrection();
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::loadRestartFiles() {
  std::cout << _systemController->getSettings().load << std::endl;
  if (_systemController->getSettings().load.empty()) {
    throw SerenityError("Option restart in SCFTask only available in combination with option load in system block!");
  }
  std::shared_ptr<OrbitalController<SCFMode>> orbitals;
  try {
    orbitals = std::make_shared<OrbitalController<SCFMode>>(_systemController->getSettings().load + "tmp",
                                                            _systemController->getBasisController(),
                                                            _systemController->getSystemIdentifier());
  }
  catch (...) {
    OutputControl::nOut << "No temporary orbital files found. Looking for converged orbital files" << std::endl;
    orbitals = std::make_shared<OrbitalController<SCFMode>>(
        _systemController->getSettings().load + _systemController->getSystemName(),
        _systemController->getBasisController(), _systemController->getSystemIdentifier());
  }
  _systemController->setElectronicStructure(std::make_shared<ElectronicStructure<SCFMode>>(
      orbitals, _systemController->getOneElectronIntegralController(), _systemController->getNOccupiedOrbitals<SCFMode>()));
}
template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::printHeader() {
  printSubSectionTitle("Main SCF Options");
  const Settings& systemSettings = _systemController->getSettings();
  auto m = systemSettings.method;
  auto s = SCFMode;
  std::string method, scfmode;
  Options::resolve<Options::ELECTRONIC_STRUCTURE_THEORIES>(method, m);
  Options::resolve<Options::SCF_MODES>(scfmode, s);
  OutputControl::n.printf("%4s SCF Mode:              %15s\n", "", scfmode.c_str());
  OutputControl::n.printf("%4s Method:                %15s\n", "", method.c_str());
  if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    if (systemSettings.customFunc.basicFunctionals.size()) {
      printf("%4s Functional:            %15s\n", "", "Custom");
      Functional(systemSettings.customFunc).print();
    }
    else {
      std::string functional;
      auto func = systemSettings.dft.functional;
      Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(functional, func);
      OutputControl::n.printf("%4s Functional:            %15s\n", "", functional.c_str());
    }
    OutputControl::n.printf("%4s Grid Accuracy:         %13d/%1d\n", "", systemSettings.grid.smallGridAccuracy,
                            systemSettings.grid.accuracy);
    std::string dispersion;
    auto disp = systemSettings.dft.dispersion;
    Options::resolve<Options::DFT_DISPERSION_CORRECTIONS>(dispersion, disp);
    OutputControl::n.printf("%4s Dispersion Correction: %15s\n", "", dispersion.c_str());
  }
  unsigned nb = _systemController->getBasisController()->getNBasisFunctions();
  OutputControl::n.printf("%4s Basis Set:             %15s\n", "", systemSettings.basis.label.c_str());
  OutputControl::n.printf("%4s Basis Functions:       %15i\n", "", nb);
  OutputControl::n.printf("%4s Caching Threshold:     %15i\n", "", systemSettings.basis.intCondition);
  std::shared_ptr<BasisController> basis = _systemController->getBasisController();
  double integralThreshold = systemSettings.basis.integralThreshold;
  if (integralThreshold == 0) {
    integralThreshold = basis->getPrescreeningThreshold();
  }
  OutputControl::n.printf("%4s Integral Threshold:    %15.1e\n", "", integralThreshold);
  OutputControl::n.printf("\n%4s Energy Threshold:      %15.1e\n", "", systemSettings.scf.energyThreshold);
  OutputControl::n.printf("%4s RMSD[D] Threshold:     %15.1e\n", "", systemSettings.scf.rmsdThreshold);
  OutputControl::n.printf("%4s DIIS Threshold:        %15.1e\n", "", systemSettings.scf.diisThreshold);
  if (_systemController->getGeometry()->hasAtomsWithECPs()) {
    OutputControl::n.printf("%4s ECP Start:             %15d\n", "", systemSettings.basis.firstECP);
  }
  if (systemSettings.pcm.use) {
    Solvents::printSolventInfo(systemSettings.pcm);
  }
  bool guess = !(_systemController->hasElectronicStructure<SCFMode>());
  if (guess) {
    auto ig = systemSettings.scf.initialguess;
    std::string init_guess;
    Options::resolve<Options::INITIAL_GUESSES>(init_guess, ig);
    OutputControl::n.printf("%4s Initial Guess:         %15s\n", "", init_guess.c_str());
  }
  std::string fittingA;
  OutputControl::n.printf("\n%4s Density Fitting:\n", "");
  auto fitj = systemSettings.basis.densFitJ;
  Options::resolve<Options::DENS_FITS>(fittingA, fitj);
  auto fitk = systemSettings.basis.densFitK;
  std::string fittingB;
  Options::resolve<Options::DENS_FITS>(fittingB, fitk);
  auto fitlrk = systemSettings.basis.densFitLRK;
  std::string fittingC;
  Options::resolve<Options::DENS_FITS>(fittingC, fitlrk);
  auto fitcorr = systemSettings.basis.densFitCorr;
  std::string fittingD;
  Options::resolve<Options::DENS_FITS>(fittingD, fitcorr);
  if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    fittingA = "NONE";
    fittingB = "NONE";
  }
  OutputControl::n.printf("%6s Coulomb              %15s\n", "", fittingA.c_str());
  OutputControl::n.printf("%6s Exchange             %15s\n", "", fittingB.c_str());
  if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    OutputControl::n.printf("%6s Long-Range Exchange  %15s\n", "", fittingC.c_str());
    OutputControl::n.printf("%6s Correlation          %15s\n", "", fittingD.c_str());
  }
  printSubSectionTitle("SCF");
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::printResults() {
  auto es = _systemController->getElectronicStructure<SCFMode>();
  auto energyComponentController = es->getEnergyComponentController();
  printSubSectionTitle("Final SCF Results");
  energyComponentController->printAllComponents();
  if (iOOptions.printFinalOrbitalEnergies) {
    print("");
    printSmallCaption("Orbital Energies");
    es->printMOEnergies();
    print("");
  }
  // SCF Analysis
  printSmallCaption("Additional Analysis");
  SCFAnalysis<SCFMode> scfAn({_systemController});
  auto s2val = scfAn.getS2();
  auto virialRatio = scfAn.getVirialRatio();
  OutputControl::n.printf("\n   -<V>/<T> = %4.3f ", virialRatio);
  OutputControl::n.printf("\n      <S*S> = %4.3f ", s2val);
  double S = fabs(0.5 * _systemController->getSpin());
  OutputControl::n.printf("\n    S*(S+1) = %4.3f ", S * (S + 1));
  OutputControl::n.printf("\n          C = %4.3f \n\n", s2val - S * (S + 1));
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED && _systemController->getSettings().scf.rohf != Options::ROHF_TYPES::NONE) {
    OutputControl::n.printf(
        "    This is a constrained UHF wavefunction without (CUHF) or reduced (SUHF)\n    spin contamination "
        "(ROHF). Take care with post HF methods or gradients.\n\n");
    if (_systemController->getSettings().scf.rohf == Options::ROHF_TYPES::SUHF) {
      OutputControl::n.printf("    Type   :     SUHF\n");
      OutputControl::n.printf("    Lambda :     %8.2e\n\n", _systemController->getSettings().scf.suhfLambda);
    }
    else if (_systemController->getSettings().scf.rohf == Options::ROHF_TYPES::CUHF) {
      OutputControl::n.printf("    Type   :     CUHF\n\n");
    }
  }
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::performSCF(std::shared_ptr<PotentialBundle<SCFMode>> potentials) {
  auto es = _systemController->getElectronicStructure<SCFMode>();
  auto energyComponentController = es->getEnergyComponentController();

  // Changing this from const reference to a copy to be able to safely turn off
  // the levelshift in case of a Delta-SCF calculation (Jul 2022, NN).
  Settings systemSettings = _systemController->getSettings();

  auto momMatrix = DeltaScf<SCFMode>::prepareMOMMatrix(settings.exca, settings.excb, es);
  // In case of a Delta-SCF calculation, turn off the levelshift
  // And make sure you're starting from a proper reference electronic structure.
  if (momMatrix != nullptr) {
    OutputControl::nOut << " Turning off levelshift for Delta-SCF calculation." << std::endl << std::endl;
    if (!es->checkFock()) {
      throw SerenityError("Make sure you have a converged ground-state electronic structure before performing a "
                          "Delta-SCF calculation.");
    }
    systemSettings.scf.useLevelshift = false;
  }
  if (this->settings.skipSCF) {
    auto F = potentials->getFockMatrix(es->getDensityMatrix(), energyComponentController);
    if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      calculateMP2Contribution();
      calculateDispersionCorrection();
    }
    es->setFockMatrix(F);
    es->attachPotentials(potentials);
    // ToDo: Crashes: EvaluateEnergyTaskTest.restricted_sDFT
    // es->toHDF5(systemSettings.path + systemSettings.name, systemSettings.identifier);
  }
  else {
    if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      // disable orbital output for dft
      auto tmp1 = iOOptions.printFinalOrbitalEnergies;
      auto tmp2 = iOOptions.gridAccuracyCheck;
      iOOptions.printFinalOrbitalEnergies = false;
      iOOptions.gridAccuracyCheck = false;

      Scf<SCFMode>::perform(systemSettings, es, potentials, settings.allowNotConverged, momMatrix, settings.momCycles);

      iOOptions.printFinalOrbitalEnergies = tmp1;
      iOOptions.gridAccuracyCheck = tmp2;

      finalDFTEnergyEvaluation(momMatrix);
    }
    else {
      Scf<SCFMode>::perform(systemSettings, es, potentials, settings.allowNotConverged, momMatrix, settings.momCycles);
    }
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<PotentialBundle<SCFMode>> ScfTask<SCFMode>::getPotentialBundle() {
  auto& systemSettings = _systemController->getSettings();
  std::shared_ptr<PotentialBundle<SCFMode>> potentials;
  if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    potentials = _systemController->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  }
  else if (systemSettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    auto gridPurpose = (this->settings.skipSCF) ? Options::GRID_PURPOSES::DEFAULT : Options::GRID_PURPOSES::SMALL;
    potentials = _systemController->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(gridPurpose);
  }
  else {
    throw SerenityError("ScfTask: Nonexistent electronicStructureTheory requested.");
  }
  return potentials;
}

template<Options::SCF_MODES SCFMode>
void ScfTask<SCFMode>::run() {
  if (this->settings.restart)
    loadRestartFiles();
  // Output of the main options
  if (iOOptions.printSCFResults)
    printHeader();
  // Get the Fock matrix routines.
  std::shared_ptr<PotentialBundle<SCFMode>> potentials = getPotentialBundle();

  // Allow fraction occupation of degenerate orbitals.
  double degeneracyThreshold = _systemController->getSettings().scf.degeneracyThreshold;
  _systemController->getElectronicStructure<SCFMode>()->getDensityMatrixController()->setDegeneracyThreshold(degeneracyThreshold);
  // Run SCF or energy evaluation.
  performSCF(potentials);

  // Print the final results.
  if (iOOptions.printSCFResults)
    printResults();
  // Clean up.
  _systemController->clear4CenterCache();
  return;
}

template class ScfTask<Options::SCF_MODES::RESTRICTED>;
template class ScfTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
