/**
 * @file   TDEmbeddingTask.cpp
 *
 * @date   Apr 23, 2014
 * @author Jan Unsleber, Thomas Dresselhaus
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
#include "tasks/TDEmbeddingTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/FockMatrix.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "grid/AtomCenteredGridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "scf/Scf.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/FDETask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemAdditionTask.h"
#include "tasks/SystemSplittingTask.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
TDEmbeddingTask<SCFMode>::TDEmbeddingTask(std::shared_ptr<SystemController> activeSystem,
                                          std::shared_ptr<SystemController> environmentSystem)
  : _activeSystem(activeSystem), _environmentSystem(environmentSystem) {
}

template<Options::SCF_MODES SCFMode>
void TDEmbeddingTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, {_activeSystem, _environmentSystem});
  // Check input for obvious misunderstanding of the keywords.
  checkInput();
  // set up the supersystem.
  std::shared_ptr<SystemController> supersystem = setUpSupersystem();
  // Get fermi-level of the supersystem.
  if (settings.useFermiLevel && settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
    settings.embedding.fermiShift = getFermiLevel(supersystem);
  }
  // Localize orbitals
  if (!settings.addOrbitals && settings.systemPartitioning != Options::SYSTEM_SPLITTING_ALGORITHM::SPADE) {
    LocalizationTask superSystemOrbLocalization(supersystem);
    superSystemOrbLocalization.settings.locType = settings.locType;
    superSystemOrbLocalization.settings.splitValenceAndCore = settings.splitValenceAndCore;
    superSystemOrbLocalization.run();
  }
  // Select orbitals
  SystemSplittingTask<SCFMode> splittingTask(supersystem, {_activeSystem, _environmentSystem});
  splittingTask.settings.systemPartitioning = settings.systemPartitioning;
  splittingTask.settings.orbitalThreshold = settings.orbitalThreshold;
  splittingTask.run();

  // Run basis-set truncation
  if (settings.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE) {
    BasisSetTruncationTask<SCFMode> basisSetTask(_activeSystem);
    basisSetTask.settings.truncationFactor = settings.truncationFactor;
    basisSetTask.settings.netThreshold = settings.netThreshold;
    basisSetTask.settings.truncAlgorithm = settings.truncAlgorithm;
    basisSetTask.run();
  }
  printSubSectionTitle("Embedded-SCF Calculation");
  if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::RECONSTRUCTION) {
    FDETask<SCFMode> fdeTask(_activeSystem, {_environmentSystem});
    fdeTask.setSuperSystemGrid(_activeSystem->getGridController());
    fdeTask.settings.embedding = settings.embedding;
    fdeTask.settings.calculateEnvironmentEnergy = true;
    fdeTask.settings.lcSettings = settings.lcSettings;
    if (settings.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE)
      fdeTask.settings.lcSettings.useProjectedOccupiedOrbitals = true;
    fdeTask.settings.maxCycles = settings.maxCycles;
    fdeTask.settings.maxResidual = settings.maxResidual;
    fdeTask.settings.loc.locType = settings.locType;
    fdeTask.settings.mp2Type = settings.mp2Type;
    fdeTask.run();
  }
  else {
    runTDPotentialReconstruction(supersystem);
  }
}

template<Options::SCF_MODES SCFMode>
void TDEmbeddingTask<SCFMode>::runTDPotentialReconstruction(std::shared_ptr<SystemController> supersystem) {
  // Turn off the system-specific continuum model.
  bool oldContinuumModelMode = _activeSystem->getSystemContinuumModelMode();
  if (oldContinuumModelMode)
    _activeSystem->setSystemContinuumModelMode(false);
  /*
   * The following code performs a FDE-type calculation.
   * The FDETask is not called here, because potential reconstruction is significantly
   * different between top-down and bottom-up-type calculations.
   */
  auto densityMatrixEnvironment = _environmentSystem->getElectronicStructure<SCFMode>()->getDensityMatrix();
  // Building the density matrix controllers needed for the potential evaluation
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMats = {
      std::make_shared<DensityMatrixController<SCFMode>>(densityMatrixEnvironment)};
  /*
   * Create Potentials
   */
  auto test = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  /*
   * Funny enough the "supersystemGrid" is built with the settings of the active system while the
   * grid of the supersystem is built with the settings of the environment.
   */
  std::shared_ptr<AtomCenteredGridController> supersystemGrid = _activeSystem->getAtomCenteredGridController();
  auto pbePot = FDEPotentialBundleFactory<SCFMode>::produce(
      _activeSystem, _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      {_environmentSystem}, envDensMats, std::make_shared<EmbeddingSettings>(settings.embedding), supersystemGrid,
      supersystem, true, settings.noSupRec);
  const auto& actsettings = _activeSystem->getSettings();
  // non-additive Dispersion correction
  auto nadDispCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      settings.embedding.dispersion, supersystem->getGeometry(), settings.embedding.naddXCFunc);
  nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      settings.embedding.dispersion, _environmentSystem->getGeometry(), settings.embedding.naddXCFunc);
  nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      settings.embedding.dispersion, _activeSystem->getGeometry(), settings.embedding.naddXCFunc);
  auto actSysDisp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      actsettings.dft.dispersion, _activeSystem->getGeometry(), _activeSystem->getSettings().dft.functional);

  /*
   * Calculate the frozen environment energy contributions and add them.
   */
  auto envEnergies = _environmentSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController();

  // add perturbative energy correction to the EnergyController.
  // TODO: If local MP2 is available, this may be different from zero.
  //       But will have to be calculate after the SCF for the active system
  //       in order to be able to project the occupied active orbitals.
  envEnergies->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, 0.0));
  std::shared_ptr<PotentialBundle<SCFMode>> envPot;
  double frozenEnvEnergy = 0.0;
  auto envsettings = _environmentSystem->getSettings();
  if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    envPot = _environmentSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
    frozenEnvEnergy = envEnergies->getEnergyComponent(ENERGY_CONTRIBUTIONS::HF_ENERGY);
  }
  else if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    // Dispersion correction
    auto envDispEnergy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
        _environmentSystem->getSettings().dft.dispersion, _environmentSystem->getGeometry(),
        _environmentSystem->getSettings().dft.functional);
    envEnergies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, envDispEnergy);
    envPot = _environmentSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
    envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
    frozenEnvEnergy = envEnergies->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
  }
  else {
    throw SerenityError("Nonexistent electronicStructureTheory requested. Options are HF and DFT.");
  }
  // Everything is set in environment ElectronicStructure -> save to file
  _environmentSystem->template getElectronicStructure<SCFMode>()->toHDF5(_environmentSystem->getHDF5BaseName(),
                                                                         _environmentSystem->getSystemIdentifier());

  auto es = _activeSystem->getElectronicStructure<SCFMode>();
  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, actSysDisp);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP, nadDispCorrection);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, frozenEnvEnergy);

  /*
   * Attach everything and run the SCF procedure
   */
  // manage output
  auto tmp1 = iOOptions.printFinalOrbitalEnergies;
  auto tmp2 = iOOptions.gridAccuracyCheck;
  iOOptions.printFinalOrbitalEnergies = false;
  iOOptions.gridAccuracyCheck = false;
  Scf<SCFMode>::perform(actsettings, es, pbePot);
  iOOptions.printFinalOrbitalEnergies = tmp1;
  iOOptions.gridAccuracyCheck = tmp2;

  // check for double hybrid functional
  auto functional = actsettings.customFunc.basicFunctionals.size() ? Functional(actsettings.customFunc)
                                                                   : resolveFunctional(actsettings.dft.functional);
  double MP2Correlation = 0.0;
  if (functional.isDoubleHybrid()) {
    WarningTracker::printWarning(
        "Warning: When using an double hybrid functional in embedding, the occupied environment orbitals are not "
        "explicitly projected out of the MP2 calculation! This may lead to errors.",
        iOOptions.printSCFCycleInfo);
    // perform MP2 for double hybrids
    if (actsettings.basis.densFitCorr != Options::DENS_FITS::NONE) {
      _activeSystem->setBasisController(
          std::dynamic_pointer_cast<AtomCenteredBasisController>(
              supersystem->getAuxBasisController(Options::AUX_BASIS_PURPOSES::CORRELATION, actsettings.basis.densFitCorr)),
          supersystem->resolveAuxBasisPurpose(Options::AUX_BASIS_PURPOSES::CORRELATION, actsettings.basis.densFitCorr));
      RIMP2<SCFMode> rimp2(_activeSystem, functional.getssScaling(), functional.getosScaling());
      MP2Correlation = rimp2.calculateCorrection();
    }
    else {
      assert(SCFMode == RESTRICTED && "MP2 is not available for unrestricted systems please use RI-MP2.");
      MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_activeSystem, functional.getssScaling(), functional.getosScaling());
      MP2Correlation = mp2EnergyCorrector.calculateElectronicEnergy();
    }
    MP2Correlation *= functional.getHfCorrelRatio();
  }

  // add energy to the EnergyController
  eCont->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, MP2Correlation));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_MP2_INT_ENERGY, 0.0);

  if (iOOptions.printSCFResults) {
    printSubSectionTitle("Final SCF Results");
    eCont->printAllComponents();
    if (iOOptions.printFinalOrbitalEnergies) {
      print("");
      printSmallCaption("Orbital Energies");
      es->printMOEnergies();
      print("");
    }
  }
  _activeSystem->setSystemContinuumModelMode(oldContinuumModelMode);
  return;
}

template<Options::SCF_MODES SCFMode>
inline void TDEmbeddingTask<SCFMode>::checkInput() {
  if (settings.embedding.carterCycles != 0 and !settings.noSupRec) {
    throw SerenityError(
        (std::string) "Zhang-Carter reconstruction does not support the double reconstruction feature!");
  }
}
template<Options::SCF_MODES SCFMode>
inline std::shared_ptr<SystemController> TDEmbeddingTask<SCFMode>::setUpSupersystem() {
  std::shared_ptr<SystemController> supersystem;
  if (settings.load == "" || settings.name == "") {
    // set up the supersystem.
    Settings supersystemSettings = _environmentSystem->getSettings();
    supersystemSettings.name = _activeSystem->getSystemName() + "+" + _environmentSystem->getSystemName();
    supersystemSettings.charge = 0;
    supersystemSettings.spin = 0;
    supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersystemSettings);
    SystemAdditionTask<SCFMode> additionTask(supersystem, {_activeSystem, _environmentSystem});
    additionTask.settings.addOccupiedOrbitals = settings.addOrbitals;
    additionTask.run();
    if (!settings.addOrbitals) {
      // run supersystem SCF
      printSubSectionTitle("Initial Supersystem-SCF Calculation");
      ScfTask<SCFMode> supersystemSCF(supersystem);
      supersystemSCF.settings.mp2Type = settings.mp2Type;
      supersystemSCF.run();
    } // if !settings.addOrbitals
  }
  else {
    OutputControl::mOut << "      Loading system: " + settings.name + " as supersystem" << std::endl;
    Settings superSysSettings;
    superSysSettings.load = settings.load;
    superSysSettings.name = settings.name;
    supersystem = std::make_shared<SystemController>(superSysSettings);
  }
  return supersystem;
}

template<Options::SCF_MODES SCFMode>
inline double TDEmbeddingTask<SCFMode>::getFermiLevel(std::shared_ptr<SystemController> supersystem) {
  std::shared_ptr<PotentialBundle<SCFMode>> potBundle;
  if (supersystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    potBundle = supersystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  }
  else if (supersystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    potBundle = supersystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  }
  else {
    throw SerenityError("Nonexistent electronicStructureTheory requested. Options are HF and DFT.");
  }
  const FockMatrix<SCFMode> fockMatrix =
      potBundle->getFockMatrix(supersystem->getElectronicStructure<SCFMode>()->getDensityMatrix(),
                               supersystem->getElectronicStructure<SCFMode>()->getEnergyComponentController());
  const CoefficientMatrix<SCFMode> coefficients = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  const auto nOcc = supersystem->getNOccupiedOrbitals<SCFMode>();
  double fermiLevel = 0.0;
  for_spin(fockMatrix, coefficients, nOcc) {
    double newMaxCoefficient =
        (coefficients_spin.leftCols(nOcc_spin).transpose() * fockMatrix_spin * coefficients_spin.leftCols(nOcc_spin))
            .diagonal()
            .array()
            .maxCoeff();
    if (newMaxCoefficient > fermiLevel)
      fermiLevel = newMaxCoefficient;
  };
  return fermiLevel;
}

template class TDEmbeddingTask<Options::SCF_MODES::RESTRICTED>;
template class TDEmbeddingTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
