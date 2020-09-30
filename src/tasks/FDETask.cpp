/**
 * @file   FDETask.cpp
 * @author Jan Unsleber, Kevin Klahr, Thomas Dresselhaus
 *
 * @date   last reworked on Nov 11. 2016
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
#include "tasks/FDETask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/MatrixInBasis.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "geometry/MolecularSurfaceFactory.h"
#include "grid/GridControllerFactory.h"
#include "io/FormattedOutput.h"
#include "misc/WarningTracker.h"
#include "postHF/MPn/LocalMP2.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "scf/Scf.h"
#include "settings/Settings.h"
#include "tasks/LocalizationTask.h"
/* Include Std and External Headers */
#include <cassert>
#include <limits>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FDETask<SCFMode>::FDETask(std::shared_ptr<SystemController> activeSystem,
                          std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems), _supersystemgrid(nullptr), _finalGrid(nullptr) {
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::run() {
  /*
   * Initialize data objects and gather data
   */
  // settings
  const auto& actsettings = _activeSystem->getSettings();

  // list of all systems for easy loops
  auto allsystems = _environmentSystems;
  allsystems.push_back(_activeSystem);

  // atoms of all subsystems
  std::vector<std::shared_ptr<Atom>> superSystemAtomsGrid;
  std::vector<std::shared_ptr<Atom>> superSystemAtomsCavity;
  superSystemAtomsCavity.insert(superSystemAtomsCavity.end(), _activeSystem->getAtoms().begin(),
                                _activeSystem->getAtoms().end());
  if (!_supersystemgrid) {
    superSystemAtomsGrid = superSystemAtomsCavity;
  }

  _activeSystem->setDiskMode(false);
  for (auto sys : _environmentSystems) {
    if (sys->getSettings().scfMode == RESTRICTED) {
      sys->template getElectronicStructure<RESTRICTED>();
    }
    else {
      sys->template getElectronicStructure<UNRESTRICTED>();
    }
    if (settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC)
      sys->setDiskMode(true);
    superSystemAtomsCavity.insert(superSystemAtomsCavity.end(), sys->getAtoms().begin(), sys->getAtoms().end());

    if (!_supersystemgrid) {
      double cutoff = settings.gridCutOff;
      if (cutoff < 0.0)
        cutoff = std::numeric_limits<double>::infinity();
      for (auto atom : sys->getAtoms()) {
        for (auto check : _activeSystem->getAtoms()) {
          double dist = distance(*atom, *check);
          if (dist < cutoff) {
            superSystemAtomsGrid.push_back(atom);
            break;
          }
        }
      }
    }
  }
  if (!_supersystemgrid) {
    // geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtomsGrid);
    superSystemGeometry->deleteIdenticalAtoms();
    // supersystem grid
    Options::GRID_PURPOSES gridacc =
        (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL : Options::GRID_PURPOSES::DEFAULT;
    _supersystemgrid = GridControllerFactory::produce(superSystemGeometry, _activeSystem->getSettings(), gridacc);
  }
  if (settings.initializeSuperMolecularSurface && settings.embedding.pcm.use) {
    // geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtomsCavity);
    superSystemGeometry->deleteIdenticalAtoms();
    auto molecularSurface = MolecularSurfaceFactory::produce(superSystemGeometry, settings.embedding.pcm);
    _activeSystem->setMolecularSurface(molecularSurface, MOLECULAR_SURFACE_TYPES::FDE);
    for (auto sys : _environmentSystems)
      sys->setMolecularSurface(molecularSurface, MOLECULAR_SURFACE_TYPES::FDE);
  }
  // Set the grid controller to the supersystem grid if "exact" methods are used.
  if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
    _activeSystem->setGridController(_supersystemgrid);
  }

  // geometries of the environment subsystems
  std::vector<std::shared_ptr<const Geometry>> envGeometries;
  for (auto sys : _environmentSystems) {
    envGeometries.push_back(sys->getGeometry());
  }

  /*
   * Initial DFT calculations on the subsystems if needed
   */
  auto es = _activeSystem->template getElectronicStructure<SCFMode>();
  /*
   * The total environment energy:
   *  - The sum of all isolated KS Energies
   *  - The additive interactions between all environment systems (ESP env-env)
   *  - Because the ESP act-env is partially inside the terms added here this part is subtracted
   *    at the start.
   */
  std::vector<std::shared_ptr<EnergyComponentController>> eConts;
  double totEnvEnergy = 0.0;
  // Check if interaction energies of active systems are present at this point
  bool allInteractionsPresent =
      _activeSystem->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY);
  eConts.push_back(_activeSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController());

  for (auto sys : _environmentSystems) {
    assert((sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
            sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) &&
           "Unknown electronic structure theory!");
    if (settings.calculateEnvironmentEnergy) {
      auto envsettings = sys->getSettings();
      auto envEnergies = sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
      auto densityMatrixEnvironment = sys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
      if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
        auto envPot = sys->template getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
        envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
      }
      else if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
        // Dispersion correction
        auto envDispEnergy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
            sys->getSettings().dft.dispersion, sys->getGeometry(), sys->getSettings().dft.functional);
        envEnergies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, envDispEnergy);
        auto envPot = sys->template getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
        envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
        envEnergies->addOrReplaceComponent(
            std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, 0.0));
        auto fBaseName = sys->getHDF5BaseName();
        if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
          fBaseName = fBaseName + ".energies.unres";
        }
        else {
          fBaseName = fBaseName + ".energies.res";
        }
        envEnergies->toFile(fBaseName, sys->getSettings().identifier);
      }
      else {
        assert(false && "Non-existing electronic structure theory requested");
      }
    } // if settings.calculateEnvironmentEnergy

    auto energyContribution = (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
                                  ? ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY
                                  : ENERGY_CONTRIBUTIONS::HF_ENERGY;
    if (sys->getSettings().scfMode == SCFMode) {
      // Add isolated energies to totEnvEnergy; we need them either way.
      if (sys->template getElectronicStructure<SCFMode>()->checkEnergy(energyContribution)) {
        totEnvEnergy += sys->template getElectronicStructure<SCFMode>()->getEnergy(energyContribution);
      }
      // Check interaction energies of every env system are present at this point
      allInteractionsPresent *=
          sys->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
          (sys->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY) &&
           sys->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
      eConts.push_back(sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    }
    else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        WarningTracker::printWarning("Warning: Calculating environment system " + sys->getSettings().name +
                                         " in restricted scfMode.\n Will proceed by assuming that alpha and beta "
                                         "densities are equal for this system...\n",
                                     iOOptions.printSCFCycleInfo);
        // Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(energyContribution)) {
          totEnvEnergy +=
              sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energyContribution);
        }
        // Check interaction energies of every env system are present at this point
        allInteractionsPresent *= sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
                                  (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(
                                       ENERGY_CONTRIBUTIONS::HF_ENERGY) &&
                                   sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(
                                       ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
        eConts.push_back(sys->template getElectronicStructure<RESTRICTED>()->getEnergyComponentController());
      }
      else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        WarningTracker::printWarning(
            "Warning: Calculating environment system " + sys->getSettings().name +
                " in unrestricted scfMode.\n Will proceed by summing up alpha and beta densities for this system...\n",
            iOOptions.printSCFCycleInfo);
        // Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(energyContribution)) {
          totEnvEnergy +=
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energyContribution);
        }
        // Check interaction energies of every env system are present at this point
        allInteractionsPresent *= sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
                                  (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(
                                       ENERGY_CONTRIBUTIONS::HF_ENERGY) &&
                                   sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(
                                       ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
        eConts.push_back(sys->template getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController());
      }
      else {
        assert(false);
      }
    }
  }
  /*
   * If all interaction energies are present:
   * Add them up (with a factor of 0.5 to prevent double counting). Then subtract
   * the interaction energy towards the active system, since it should be part
   * of the active systems energy (also with a factor 0.5, since we added only
   * half it in the first place...
   */
  if (allInteractionsPresent) {
    for (auto sys : _environmentSystems) {
      if (sys->getSettings().scfMode == SCFMode) {
        totEnvEnergy += 0.5 * sys->template getElectronicStructure<SCFMode>()->getEnergy(
                                  ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
      }
      else {
        if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
        }
        else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
        }
        else {
          assert(false);
        }
      }
    }
    totEnvEnergy -= 0.5 * _activeSystem->template getElectronicStructure<SCFMode>()->getEnergy(
                              ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
  }
  else if (!settings.calculateEnvironmentEnergy) {
    WarningTracker::printWarning("WARNING: Not all interaction energies available. The environment energy is\n \
        calculated using only isolated KS-DFT energies. Be aware that this does not lead to\n \
        the correct supersystem energy. If you want the supersystem energy, please calculate\n \
        the interaction energies for every subsystem (e.g., by choosing them as active system\n \
        in a FDE or FAT calculation.)\n",
                                 iOOptions.printSCFCycleInfo);
  }

  // list of environment density matrices (their controllers)
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensities;
  for (auto sys : _environmentSystems) {
    if (sys->getSettings().scfMode == SCFMode) {
      envDensities.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    }
    else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        // Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      }
      else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        // Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin =
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      }
      else {
        assert(false);
      }
    }
  }

  // Turn off the system-specific continuum model.
  bool oldContinuumModelMode = _activeSystem->getSystemContinuumModelMode();
  if (oldContinuumModelMode)
    _activeSystem->setSystemContinuumModelMode(false);

  /*
   * Create Potentials
   */
  auto fdePot = FDEPotentialBundleFactory<SCFMode>::produce(
      _activeSystem, _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      _environmentSystems, envDensities, std::make_shared<EmbeddingSettings>(settings.embedding), _supersystemgrid,
      nullptr, false, true, settings.gridCutOff, eConts, settings.firstPassiveSystemIndex);
  /*
   * non-additive dispersion correction
   */
  calculateNonAdditiveDispersionCorrection();
  double actDispCorrection = 0.0;
  if (actsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    // Active dispersion correction.
    actDispCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
        actsettings.dft.dispersion, _activeSystem->getGeometry(), _activeSystem->getSettings().dft.functional);
  }
  /*
   * Attach everything and run scf
   */
  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, actDispCorrection);
  /*
   * Attach everything and run scf
   */
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, totEnvEnergy);
  Scf<SCFMode>::perform(actsettings, es, fdePot);

  if (settings.gridCutOff > 0.0 and settings.finalGrid and !_finalGrid) {
    // atoms of all subsystems
    auto finalGridGeometry = std::make_shared<Geometry>();
    *finalGridGeometry += *_activeSystem->getGeometry();

    for (auto sys : _environmentSystems) {
      *finalGridGeometry += *sys->getGeometry();
    }

    // supersystem grid
    Options::GRID_PURPOSES finalGridacc =
        (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL : Options::GRID_PURPOSES::DEFAULT;
    _finalGrid = GridControllerFactory::produce(finalGridGeometry, _activeSystem->getSettings(), finalGridacc);
  }

  if (settings.calculateSolvationEnergy) {
    const double actElEnvEl = eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_ELECTRONS_COULOMB);
    const double actElEnvNuc = eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_NUCLEI_COULOMB);
    const double actNucEnvEl = eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_NUCLEI_ENV_ELECTRONS_COULOMB);
    const double actNucEnvNuc = eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_NUCLEI_ENV_NUCLEI_COULOMB);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_ELECTRONS_ENV_ELECTRONS_COULOMB, 0.5 * actElEnvEl);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_ELECTRONS_ENV_NUCLEI_COULOMB, 0.5 * actElEnvNuc);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NUCLEI_ENV_ELECTRONS_COULOMB, 0.5 * actNucEnvEl);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NUCLEI_ENV_NUCLEI_COULOMB, 0.5 * actNucEnvNuc);
  }
  if (settings.gridCutOff > 0.0 and settings.finalGrid) {
    auto naddXCfunc = resolveFunctional(settings.embedding.naddXCFunc);
    auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid, naddXCfunc, std::make_pair(true, eConts), true, true);
    auto kin = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(), envDensities,
        _finalGrid, resolveFunctional(settings.embedding.naddKinFunc), std::make_pair(false, eConts), true, true);
    const auto& P =
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    double eNaddKin = kin->getEnergy(P);
    double eNaddXC = naddXCPot->getEnergy(P);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC, eNaddXC);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC, eNaddKin);
  }

  if (settings.calculateSolvationEnergy) {
    if (!_finalGrid)
      _finalGrid = _supersystemgrid;
    auto naddXCfunc = resolveFunctional(settings.embedding.naddXCFunc);
    auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid, naddXCfunc, std::make_pair(true, eConts), true, true, settings.calculateSolvationEnergy);
    auto kin = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid, resolveFunctional(settings.embedding.naddKinFunc), std::make_pair(false, eConts),
        true, true, settings.calculateSolvationEnergy);

    const auto& P =
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    double eLinNaddKin = 0.0;
    double eLinNaddXC = 0.0;
    eLinNaddKin = kin->getLinearizedEnergy(P, 0.6);
    eLinNaddXC = naddXCPot->getLinearizedEnergy(P, 0.75);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_KINETIC, eLinNaddKin);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_XC, eLinNaddXC);
  }

  auto functional = resolveFunctional(_activeSystem->getSettings().dft.functional);
  double MP2Correlation = 0.0;
  if (settings.calculateMP2Energy && functional.isDoubleHybrid()) {
    switch (settings.mp2Type) {
      case Options::MP2_TYPES::LOCAL: {
        LocalizationTask locTask(_activeSystem);
        locTask.settings.locType = settings.locType;
        locTask.settings.splitValenceAndCore = true;
        locTask.run();
        // perform MP2 for double hybrids
        settings.lcSettings.embeddingSettings = settings.embedding;
        FockMatrix<SCFMode> f = fdePot->getFockMatrix(_activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix(),
                                                      std::make_shared<EnergyComponentController>());
        auto fock = std::make_shared<FockMatrix<RESTRICTED>>(_activeSystem->getBasisController());
        for_spin(f) {
          *fock = f_spin;
        };
        auto localCorrelationController =
            std::make_shared<LocalCorrelationController>(_activeSystem, settings.lcSettings, _environmentSystems, fock);
        LocalMP2 localMP2(localCorrelationController);
        localMP2.settings.maxResidual = settings.maxResidual;
        localMP2.settings.maxCycles = settings.maxCycles;
        // ss and os scaling.
        localMP2.settings.osScaling = functional.getosScaling();
        localMP2.settings.ssScaling = functional.getssScaling();
        MP2Correlation = functional.getHfCorrelRatio() * localMP2.calculateEnergyCorrection().sum();
        break;
      }
      case Options::MP2_TYPES::RI: {
        RIMP2<SCFMode> rimp2(_activeSystem, functional.getssScaling(), functional.getosScaling());
        MP2Correlation = functional.getHfCorrelRatio() * rimp2.calculateCorrection();
        break;
      }
      case Options::MP2_TYPES::AO: {
        if (SCFMode != RESTRICTED)
          throw SerenityError("MP2 is not available for unrestricted systems. Please use the RI approximation.");
        MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_activeSystem, functional.getssScaling(), functional.getosScaling());
        MP2Correlation = functional.getHfCorrelRatio() * mp2EnergyCorrector.calculateElectronicEnergy();
        break;
      }
    }
  } // if settings.calculateMP2Energy && functional.isDoubleHybrid()
  eCont->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, MP2Correlation));
  /*
   * Print final results if sensible
   */
  auto fBaseName = _activeSystem->getHDF5BaseName();
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
    fBaseName = fBaseName + ".energies.unres";
  }
  else {
    fBaseName = fBaseName + ".energies.res";
  }
  eCont->toFile(fBaseName, _activeSystem->getSettings().identifier);
  if ((iOOptions.printSCFResults and not(settings.gridCutOff > 0.0 and !settings.finalGrid)) ||
      settings.calculateSolvationEnergy) {
    printSubSectionTitle("Final SCF Results");
    eCont->printAllComponents();
    if (iOOptions.printFinalOrbitalEnergies) {
      print("");
      printSmallCaption("Orbital Energies");
      es->printMOEnergies();
      print("");
    }
  }
  // Resett continuum model mode
  _activeSystem->setSystemContinuumModelMode(oldContinuumModelMode);
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::calculateNonAdditiveDispersionCorrection() {
  double nadDispCorrection = 0.0;
  double solvScaledNadDispCorrection = 0.0;
  if (settings.embedding.dispersion != Options::DFT_DISPERSION_CORRECTIONS::NONE) {
    Timings::takeTime("FDE -    Non-Add. Disper.");
    auto environmentGeometry = std::make_shared<Geometry>();
    for (auto sys : _environmentSystems)
      *environmentGeometry += *sys->getGeometry();
    environmentGeometry->deleteIdenticalAtoms();
    if (settings.calculateSolvationEnergy) {
      double scaling = 0.75;
      solvScaledNadDispCorrection = scaling * DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
                                                  settings.embedding.dispersion, _activeSystem->getGeometry(),
                                                  environmentGeometry, settings.embedding.naddXCFunc);
    }
    if (settings.finalGrid) {
      auto supersystemGeometry = std::make_shared<Geometry>();
      *supersystemGeometry += *_activeSystem->getGeometry();
      *supersystemGeometry += *environmentGeometry;
      supersystemGeometry->deleteIdenticalAtoms();
      nadDispCorrection += DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
          settings.embedding.dispersion, supersystemGeometry, settings.embedding.naddXCFunc);
      for (auto sys : _environmentSystems)
        nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
            settings.embedding.dispersion, sys->getGeometry(), settings.embedding.naddXCFunc);
    } // if settings.finalGrid
    Timings::timeTaken("FDE -    Non-Add. Disper.");
  }
  auto eCont = _activeSystem->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  if (settings.calculateSolvationEnergy) {
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_DISP, solvScaledNadDispCorrection);
  }
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP, nadDispCorrection);
}

template class FDETask<Options::SCF_MODES::RESTRICTED>;
template class FDETask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
