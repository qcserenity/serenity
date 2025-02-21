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
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "geometry/MolecularSurfaceController.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h"
#include "postHF/MPn/LTMP2.h"
#include "postHF/MPn/LocalMP2.h"
#include "postHF/MPn/LocalMP2InteractionCalculator.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "potentials/ExchangeInteractionPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "scf/SCFAnalysis.h"
#include "scf/Scf.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <limits>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FDETask<SCFMode>::FDETask(std::shared_ptr<SystemController> activeSystem,
                          std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems), _supersystemgrid(nullptr), _finalGrid(nullptr) {
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _activeSystem);
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
    if (sys->getSCFMode() == RESTRICTED) {
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
    _supersystemgrid =
        AtomCenteredGridControllerFactory::produce(superSystemGeometry, _activeSystem->getSettings().grid, gridacc);
  }
  if (settings.initializeSuperMolecularSurface && settings.embedding.pcm.use && !settings.embedding.pcm.loadedPCM) {
    // geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtomsCavity);
    superSystemGeometry->deleteIdenticalAtoms();
    auto molecularSurface = std::make_shared<MolecularSurfaceController>(superSystemGeometry, settings.embedding.pcm);
    _activeSystem->setMolecularSurface(molecularSurface, MOLECULAR_SURFACE_TYPES::FDE);
    for (auto sys : _environmentSystems)
      sys->setMolecularSurface(molecularSurface, MOLECULAR_SURFACE_TYPES::FDE);
  }
  else if (settings.embedding.pcm.use && settings.embedding.pcm.loadedPCM) {
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtomsCavity);
    superSystemGeometry->deleteIdenticalAtoms();
    auto molecularSurface = _activeSystem->getMolecularSurface(MOLECULAR_SURFACE_TYPES::FDE);
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
      _activeSystem->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
      _activeSystem->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_HF_ENERGY);
  eConts.push_back(_activeSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController());

  for (auto sys : _environmentSystems) {
    if (settings.calculateEnvironmentEnergy) {
      if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        auto envsettings = sys->getSettings();
        auto envEnergies =
            sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController();
        auto densityMatrixEnvironment =
            sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
          auto envPot =
              sys->template getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
          envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
        }
        else if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
          // Dispersion correction
          auto envDispEnergy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
              sys->getSettings().dft.dispersion, sys->getGeometry(), sys->getSettings().dft.functional);
          envEnergies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, envDispEnergy);
          auto envPot =
              sys->template getPotentials<Options::SCF_MODES::RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
          envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
          envEnergies->addOrReplaceComponent(
              std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, 0.0));
          auto fBaseName = sys->getHDF5BaseName() + ".energies.res";
          envEnergies->toFile(fBaseName, sys->getSystemIdentifier());
        }
      } // if scfmode
      else {
        auto envsettings = sys->getSettings();
        auto envEnergies =
            sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergyComponentController();
        auto densityMatrixEnvironment =
            sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix();
        if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
          auto envPot =
              sys->template getPotentials<Options::SCF_MODES::UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
          envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
        }
        else if (envsettings.method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
          // Dispersion correction
          auto envDispEnergy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
              sys->getSettings().dft.dispersion, sys->getGeometry(), sys->getSettings().dft.functional);
          envEnergies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, envDispEnergy);
          auto envPot =
              sys->template getPotentials<Options::SCF_MODES::UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
          envPot->getFockMatrix(densityMatrixEnvironment, envEnergies);
          envEnergies->addOrReplaceComponent(
              std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, 0.0));
          auto fBaseName = sys->getHDF5BaseName() + ".energies.unres";
          envEnergies->toFile(fBaseName, sys->getSystemIdentifier());
        }
      } // else scfmode
    }   // if settings.calculateEnvironmentEnergy

    auto energyContribution = (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
                                  ? ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY
                                  : ENERGY_CONTRIBUTIONS::HF_ENERGY;
    if (sys->getSCFMode() == SCFMode) {
      // Add isolated energies to totEnvEnergy; we need them either way.
      if (sys->template getElectronicStructure<SCFMode>()->checkEnergy(energyContribution)) {
        totEnvEnergy += sys->template getElectronicStructure<SCFMode>()->getEnergy(energyContribution);
      }
      // Check interaction energies of every env system are present at this point
      allInteractionsPresent =
          allInteractionsPresent &&
          (sys->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
           (sys->template getElectronicStructure<SCFMode>()->checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_HF_ENERGY)));
      eConts.push_back(sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    }
    else {
      if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        WarningTracker::printWarning("Warning: Calculating environment system " + sys->getSystemName() +
                                         " in restricted scfMode.\n Will proceed by assuming that alpha and beta "
                                         "densities are equal for this system...\n",
                                     iOOptions.printSCFCycleInfo);
        // Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(energyContribution)) {
          totEnvEnergy +=
              sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energyContribution);
        }
        // Check interaction energies of every env system are present at this point
        allInteractionsPresent = allInteractionsPresent &&
                                 (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
                                  (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_HF_ENERGY)));
        eConts.push_back(sys->template getElectronicStructure<RESTRICTED>()->getEnergyComponentController());
      }
      else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        WarningTracker::printWarning(
            "Warning: Calculating environment system " + sys->getSystemName() +
                " in unrestricted scfMode.\n Will proceed by summing up alpha and beta densities for this system...\n",
            iOOptions.printSCFCycleInfo);
        // Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(energyContribution)) {
          totEnvEnergy +=
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energyContribution);
        }
        // Check interaction energies of every env system are present at this point
        allInteractionsPresent = allInteractionsPresent &&
                                 (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY) ||
                                  (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(
                                      ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_HF_ENERGY)));
        eConts.push_back(sys->template getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController());
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
      if (sys->getSCFMode() == SCFMode) {
        totEnvEnergy += 0.5 * sys->template getElectronicStructure<SCFMode>()->getEnergy(
                                  ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
      }
      else {
        if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
        }
        else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(
                                    ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
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
    if (sys->getSCFMode() == SCFMode) {
      envDensities.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    }
    else {
      if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        // Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      }
      else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        // Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin =
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
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
  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, actDispCorrection);
  /*
   * Initiate ALMO procedure
   */
  if (settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::ALMO) {
    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;
    MatrixInBasis<SCFMode> customS(_activeSystem->getBasisController());
    Eigen::MatrixXd S =
        _activeSystem->getElectronicStructure<SCFMode>()->getOneElectronIntegralController()->getOverlapIntegrals();
    for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
      auto& libint = Libint::getInstance();
      Eigen::MatrixXd S_AB = Eigen::MatrixXd(libint.compute1eInts(
          LIBINT_OPERATOR::overlap, _environmentSystems[iEnv]->getBasisController(), _activeSystem->getBasisController()));
      const DensityMatrix<SCFMode> D_BB = _environmentSystems[iEnv]->getElectronicStructure<SCFMode>()->getDensityMatrix();
      for_spin(customS, D_BB) {
        customS_spin = S - scfFactor * S_AB * D_BB_spin * S_AB.transpose();
      };
    }
    _activeSystem->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->useCustomOverlap(customS, true);
  }
  /*
   * Attach everything and run scf
   */
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, totEnvEnergy);
  if (this->settings.skipSCF) {
    auto F = fdePot->getFockMatrix(es->getDensityMatrix(), eCont);
    es->attachPotentials(fdePot);
    _activeSystem->template getElectronicStructure<SCFMode>()->setFockMatrix(F);
  }
  else {
    Scf<SCFMode>::perform(actsettings, es, fdePot, false, nullptr, 0,
                          settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::ALMO);
  }

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
    _finalGrid =
        AtomCenteredGridControllerFactory::produce(finalGridGeometry, _activeSystem->getSettings().grid, finalGridacc);
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
    auto naddXCfunc = settings.embedding.customNaddXCFunc.basicFunctionals.size()
                          ? Functional(settings.embedding.customNaddXCFunc)
                          : resolveFunctional(settings.embedding.naddXCFunc);
    auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid, naddXCfunc, std::make_pair(true, eConts), true, true);
    auto kin = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid,
        settings.embedding.customNaddKinFunc.basicFunctionals.size() ? Functional(settings.embedding.customNaddKinFunc)
                                                                     : resolveFunctional(settings.embedding.naddKinFunc),
        std::make_pair(false, eConts), true, true);
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
    auto naddXCfunc = settings.embedding.customNaddXCFunc.basicFunctionals.size()
                          ? Functional(settings.embedding.customNaddXCFunc)
                          : resolveFunctional(settings.embedding.naddXCFunc);
    auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid, naddXCfunc, std::make_pair(true, eConts), true, true, settings.calculateSolvationEnergy);
    auto kin = std::make_shared<NAddFuncPotential<SCFMode>>(
        _activeSystem, _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities, _finalGrid,
        settings.embedding.customNaddKinFunc.basicFunctionals.size() ? Functional(settings.embedding.customNaddKinFunc)
                                                                     : resolveFunctional(settings.embedding.naddKinFunc),
        std::make_pair(false, eConts), true, true, settings.calculateSolvationEnergy);

    const auto& P =
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    double eLinNaddKin = 0.0;
    double eLinNaddXC = 0.0;
    eLinNaddKin = kin->getLinearizedEnergy(P, 0.6);
    eLinNaddXC = naddXCPot->getLinearizedEnergy(P, 0.75);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_KINETIC, eLinNaddKin);
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_XC, eLinNaddXC);
  }
  calculateMP2CorrelationContribution(fdePot);
  if (settings.calculateUnrelaxedMP2Density)
    calculateUnrelaxedMP2Density(fdePot, eCont);
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
  eCont->toFile(fBaseName, _activeSystem->getSystemIdentifier());
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
    this->printFaTAnalysis({_activeSystem}, nullptr, true);

    std::vector<std::shared_ptr<SystemController>> Systems;
    Systems.push_back(_activeSystem);
    for (auto& sys : _environmentSystems) {
      Systems.push_back(sys);
    }

    if (!_finalGrid) {
      this->printFaTAnalysis(Systems, _supersystemgrid);
    }
    else {
      this->printFaTAnalysis(Systems, _finalGrid);
    }
  }
  // Reset continuum model mode
  _activeSystem->setSystemContinuumModelMode(oldContinuumModelMode);
  // Clean up.
  _activeSystem->clear4CenterCache();
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::calculateMP2CorrelationContribution(std::shared_ptr<PotentialBundle<SCFMode>> fdePot) {
  auto eCont = _activeSystem->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  auto functional = _activeSystem->getSettings().customFunc.basicFunctionals.size()
                        ? Functional(_activeSystem->getSettings().customFunc)
                        : resolveFunctional(_activeSystem->getSettings().dft.functional);
  auto naddXCFunctional = settings.embedding.customNaddXCFunc.basicFunctionals.size()
                              ? Functional(settings.embedding.customNaddXCFunc)
                              : resolveFunctional(settings.embedding.naddXCFunc);
  double MP2Correlation = 0.0;
  double mp2Interaction = 0.0;
  if (settings.calculateMP2Energy) {
    if (functional.isDoubleHybrid() and _activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      switch (settings.mp2Type) {
        case Options::MP2_TYPES::LOCAL: {
          LocalizationTask locTask(_activeSystem);
          locTask.settings = settings.loc;
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
        case Options::MP2_TYPES::DF: {
          if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC)
            WarningTracker::printWarning(
                "WARNING: The virtual orbital space contains occupied environment orbitals.\n"
                "         This may lead to an error. Use \"mp2Type local\" to avoid this issue.",
                iOOptions.printSCFCycleInfo);
          RIMP2<SCFMode> rimp2(_activeSystem, functional.getssScaling(), functional.getosScaling());
          MP2Correlation = functional.getHfCorrelRatio() * rimp2.calculateCorrection();
          break;
        }
        case Options::MP2_TYPES::LT: {
          LTMP2<SCFMode> ltmp2(_activeSystem);
          MP2Correlation = functional.getHfCorrelRatio() * ltmp2.calculateCorrection();
          break;
        }
        case Options::MP2_TYPES::AO: {
          if (SCFMode != RESTRICTED)
            throw SerenityError("MP2 is not available for unrestricted systems. Please use the RI approximation.");
          MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_activeSystem, functional.getssScaling(),
                                                            functional.getosScaling());
          MP2Correlation = functional.getHfCorrelRatio() * mp2EnergyCorrector.calculateElectronicEnergy();
          break;
        }
      }
    } // if functional.isDoubleHybrid()
    if (naddXCFunctional.isDoubleHybrid()) {
      LocalizationTask locTaskAct(_activeSystem);
      locTaskAct.settings = settings.loc;
      locTaskAct.run();
      for (auto sys : _environmentSystems) {
        LocalizationTask locTaskEnv(sys);
        locTaskEnv.settings = settings.loc;
        locTaskEnv.run();
      }
      settings.lcSettings.embeddingSettings = settings.embedding;
      FockMatrix<SCFMode> f = fdePot->getFockMatrix(_activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix(),
                                                    std::make_shared<EnergyComponentController>());
      auto fock = std::make_shared<FockMatrix<RESTRICTED>>(_activeSystem->getBasisController());
      for_spin(f) {
        *fock = f_spin;
      };
      LocalMP2InteractionCalculator locMP2Inter(_activeSystem, _environmentSystems, settings.lcSettings, settings.maxCycles,
                                                settings.maxResidual, settings.embedding.fullMP2Coupling);
      auto envXCFunctional = _environmentSystems[0]->getSettings().customFunc.basicFunctionals.size()
                                 ? Functional(_environmentSystems[0]->getSettings().customFunc)
                                 : resolveFunctional(_environmentSystems[0]->getSettings().dft.functional);
      bool envIsDFT = _environmentSystems[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
      if (settings.calculateEnvironmentEnergy && envXCFunctional.isDoubleHybrid() && envIsDFT) {
        double envEnergy = eCont->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES);
        double envMP2Energy = locMP2Inter.getEnvironmentEnergy();
        envEnergy += envMP2Energy;
        auto envECont = _environmentSystems[0]->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
        eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES, envEnergy);
        envECont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, envMP2Energy);
      }
      mp2Interaction = locMP2Inter.getInteractionEnergy();
      mp2Interaction += locMP2Inter.getCouplingEnergyCorrection();
    } // if naddXCFunctional.isDoubleHybrid()
  }   // if settings.calculateMP2Energy
  eCont->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS, double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION, MP2Correlation));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_MP2_INT_ENERGY, mp2Interaction);
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
      nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
          settings.embedding.dispersion, _activeSystem->getGeometry(), settings.embedding.naddXCFunc);
    } // if settings.finalGrid
    Timings::timeTaken("FDE -    Non-Add. Disper.");
  }
  auto eCont = _activeSystem->getElectronicStructure<SCFMode>()->getEnergyComponentController();
  if (settings.calculateSolvationEnergy) {
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_DISP, solvScaledNadDispCorrection);
  }
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP, nadDispCorrection);
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::calculateUnrelaxedMP2Density(std::shared_ptr<PotentialBundle<SCFMode>> fdePot,
                                                    std::shared_ptr<EnergyComponentController> eCont) {
  DensityMatrix<SCFMode> densityCorrection(_activeSystem->getBasisController());
  double MP2Correlation = 0.0;
  switch (settings.mp2Type) {
    case Options::MP2_TYPES::LOCAL: {
      if (SCFMode != RESTRICTED)
        throw SerenityError("DLPNO-MP2 is not available for unrestricted systems. Please use the RI approximation.");
      LocalizationTask locTask(_activeSystem);
      locTask.settings = settings.loc;
      locTask.run();
      auto localCorrelationController =
          std::make_shared<LocalCorrelationController>(_activeSystem, settings.lcSettings, _environmentSystems);
      LocalMP2 localMP2(localCorrelationController);
      localMP2.settings.maxResidual = settings.maxResidual;
      localMP2.settings.maxCycles = settings.maxCycles;

      MP2Correlation = localMP2.calculateEnergyCorrection().sum();
      densityCorrection = localMP2.calculateDensityCorrection();
      break;
    }
    case Options::MP2_TYPES::DF: {
      if (settings.embedding.embeddingMode != Options::KIN_EMBEDDING_MODES::NADD_FUNC)
        WarningTracker::printWarning("WARNING: The virtual orbital space contains occupied environment orbitals.\n"
                                     "         This may lead to an error. Use \"mp2Type local\" to avoid this issue.",
                                     iOOptions.printSCFCycleInfo);
      RIMP2<SCFMode> rimp2(_activeSystem);
      MP2Correlation = rimp2.calculateCorrection();
      densityCorrection = rimp2.calculateDensityCorrection();
      break;
    }
    case Options::MP2_TYPES::LT: {
      throw SerenityError(
          "Unrelaxed density is not available for Laplace-Transform MP2. Please use the RI or DLPNO approximation.");
    }
    case Options::MP2_TYPES::AO: {
      throw SerenityError(
          "Unrelaxed density is not available for non-RI MP2. Please use the RI or DLPNO approximation.");
    }
  }
  auto densityController = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  auto density = densityController->getDensityMatrix();
  densityController->setDensityMatrix(density + densityCorrection);
  auto eCont_tmp = std::make_shared<EnergyComponentController>();
  fdePot->getFockMatrix(densityController->getDensityMatrix(), eCont_tmp);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC,
                               eCont_tmp->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC,
                               eCont_tmp->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_ELECTRONS_COULOMB,
                               eCont_tmp->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_ELECTRONS_COULOMB));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_NUCLEI_COULOMB,
                               eCont_tmp->getEnergyComponent(ENERGY_CONTRIBUTIONS::FDE_ELECTRONS_ENV_NUCLEI_COULOMB));
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::MP2_CORRECTION, MP2Correlation);
}

template<Options::SCF_MODES SCFMode>
void FDETask<SCFMode>::printFaTAnalysis(std::vector<std::shared_ptr<SystemController>> systems,
                                        std::shared_ptr<GridController> supersystemGrid, bool actOnly) {
  bool skipBecauseOfMissingImplementation =
      !actOnly && systems[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF && SCFMode == UNRESTRICTED;
  if (skipBecauseOfMissingImplementation) {
    printf("\n");
    return;
  }
  SCFAnalysis<SCFMode> FaTanalysisSuperSys(systems, supersystemGrid);
  auto S2 = FaTanalysisSuperSys.getS2();
  if (actOnly) {
    printSmallCaption("Analysis:");
    OutputControl::n.printf("ActiveSystem <S*S> = %4.3f \n", S2);
  }
  else {
    OutputControl::n.printf("Supersystem <S*S> = %4.3f \n\n", S2);
  }
}

template class FDETask<Options::SCF_MODES::RESTRICTED>;
template class FDETask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
