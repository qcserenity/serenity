/**
 * @file FDEPotentialBundleFactory.cpp
 *
 * @author Moritz Bensberg
 * @date Aug 21, 2019
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
#include "potentials/bundles/FDEPotentialBundleFactory.h"

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "dft/Functional.h"
#include "geometry/Geometry.h"
#include "geometry/MolecularSurfaceController.h"
#include "grid/GridController.h"
#include "io/FormattedOutputStream.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

/* Fock matrix construction */
#include "potentials/BUReconstructionPotential.h"
#include "potentials/ECPInteractionPotential.h"
#include "potentials/HoffmannProjectionPotential.h"
#include "potentials/HuzinagaFDEProjectionPotential.h"
#include "potentials/LevelshiftHybridPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/PCMPotential.h"
#include "potentials/TDReconstructionPotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/bundles/PBEPotentials.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::shared_ptr<PotentialBundle<SCFMode>> FDEPotentialBundleFactory<SCFMode>::produce(
    std::shared_ptr<SystemController> activeSystem, std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
    std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex) {
  return produceNew(activeSystem, activeDensMatController, environmentSystems, envDensMatController, settings, grid,
                    supersystem, topDown, noSuperRecon, gridCutOff, eConts, firstPassiveSystemIndex);
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<PotentialBundle<SCFMode>> FDEPotentialBundleFactory<SCFMode>::produceNew(
    std::shared_ptr<SystemController> activeSystem, std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
    std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex) {
  std::shared_ptr<PotentialBundle<SCFMode>> potBundle;
  if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    potBundle = activeSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  }
  else if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    potBundle = activeSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  }
  else {
    throw SerenityError("None existing electronicStructureTheory requested. Options are HF and DFT.");
  }
  // The collected atoms and geometries will be needed further down for ECPs and Coulomb interaction.
  std::vector<std::shared_ptr<const Geometry>> environmentGeometries;
  std::vector<std::shared_ptr<BasisController>> environmentAuxilliaryBasisSets;
  std::vector<std::shared_ptr<Atom>> environmentAtoms;
  for (auto env : environmentSystems) {
    environmentGeometries.push_back(env->getGeometry());
    environmentAuxilliaryBasisSets.push_back(env->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
  }
  for (auto envGeom : environmentGeometries)
    environmentAtoms.insert(environmentAtoms.begin(), envGeom->getAtoms().begin(), envGeom->getAtoms().end());
  std::shared_ptr<PotentialBundle<SCFMode>> esiPot;
  if (activeSystem->getSettings().dft.densityFitting == Options::DENS_FITS::RI) {
    esiPot = std::make_shared<ESIPotentials<SCFMode>>(
        activeSystem, environmentSystems, activeDensMatController, activeSystem->getGeometry(), envDensMatController,
        environmentGeometries, activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
        environmentAuxilliaryBasisSets, firstPassiveSystemIndex);
  }
  else {
    esiPot = std::make_shared<ESIPotentials<SCFMode>>(activeSystem, environmentSystems, activeDensMatController,
                                                      activeSystem->getGeometry(), envDensMatController,
                                                      environmentGeometries, firstPassiveSystemIndex);
  }
  auto naddXCfunc = resolveFunctional(settings->naddXCFunc);
  // Issue a warning if exact exchange should be evaluated with non-orthogonal orbitals
  // without any correction respecting this.
  if (naddXCfunc.isHybrid() and (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC ||
                                 settings->embeddingMode == Options::KIN_EMBEDDING_MODES::NONE)) {
    WarningTracker::printWarning("Warning: Exact exchange is calculated with orbitals that may not be orthogonal.",
                                 iOOptions.printSCFCycleInfo);
  }
  auto nonAddExchange = std::make_shared<NAddFuncPotential<SCFMode>>(
      activeSystem, activeDensMatController, envDensMatController, grid, naddXCfunc,
      std::make_pair(true, (gridCutOff < 0) ? eConts : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
      (topDown || gridCutOff < 0.0) ? true : false,
      (settings->embeddingMode != Options::KIN_EMBEDDING_MODES::RECONSTRUCTION && gridCutOff < 0.0) ? true : false);
  auto ecpInt = std::make_shared<ECPInteractionPotential<SCFMode>>(activeSystem, activeSystem->getGeometry()->getAtoms(),
                                                                   environmentAtoms, envDensMatController,
                                                                   activeSystem->getBasisController());
  bool usesPCM = settings->pcm.use;
  if (activeSystem->getSettings().pcm.use) {
    WarningTracker::printWarning("Warning: An implicit solvent model is used in the active system in combination with "
                                 "explicit environment molecules! The molecular cavity may not enclose all molecules "
                                 "and the electric potential may not be calculated correctly.",
                                 iOOptions.printSCFCycleInfo);
  }
  std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> envElecPots = {};
  if (usesPCM) {
    for (auto envSys : environmentSystems) {
      envElecPots.push_back(
          envSys->getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::FDE));
    }
  }
  auto pcm = std::make_shared<PCMPotential<SCFMode>>(
      settings->pcm, activeSystem->getBasisController(), activeSystem->getGeometry(),
      (usesPCM) ? activeSystem->getMolecularSurface(MOLECULAR_SURFACE_TYPES::FDE) : nullptr,
      (usesPCM) ? activeSystem->getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::FDE)
                : nullptr,
      envElecPots);
  std::shared_ptr<Potential<SCFMode>> kin;
  if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
    kin = std::make_shared<NAddFuncPotential<SCFMode>>(
        activeSystem, activeDensMatController, envDensMatController, grid, resolveFunctional(settings->naddKinFunc),
        std::make_pair(false, (gridCutOff < 0) ? eConts : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
        (gridCutOff < 0.0) ? true : false);
  }
  else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
    kin = std::make_shared<LevelshiftHybridPotential<SCFMode>>(
        activeSystem, environmentSystems, settings->levelShiftParameter,
        (settings->basisFunctionRatio == 0.0) ? false : true, settings->basisFunctionRatio,
        settings->borderAtomThreshold, grid, resolveFunctional(settings->longRangeNaddKinFunc), topDown,
        std::make_pair(false, std::vector<std::shared_ptr<EnergyComponentController>>(0)));
  }
  else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
    auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(new FDEPotentials<SCFMode>(
        potBundle, esiPot, nonAddExchange,
        std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController())), ecpInt, pcm));
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings,
                                                                    superPot_AA, topDown, grid, gridCutOff, eConts);
  }
  else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
    auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(new FDEPotentials<SCFMode>(
        potBundle, esiPot, nonAddExchange,
        std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController())), ecpInt, pcm));
    kin = std::make_shared<HoffmannProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings,
                                                                 superPot_AA, topDown, grid, gridCutOff, eConts);
  }
  else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION) {
    if (topDown) {
      kin = std::make_shared<TDReconstructionPotential<SCFMode>>(
          activeSystem, supersystem, environmentSystems, settings->smoothFactor,
          settings->potentialBasis.empty() ? environmentSystems[0]->getBasisController()->getBasisString() : settings->potentialBasis,
          settings->singValThreshold, settings->lbDamping, settings->lbCycles, settings->carterCycles, noSuperRecon);
    }
    else {
      kin = std::make_shared<BUReconstructionPotential<SCFMode>>(
          activeSystem, activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
          envDensMatController, grid, naddXCfunc, environmentSystems, settings->smoothFactor,
          settings->potentialBasis.empty() ? activeSystem->getBasisController()->getBasisString() : settings->potentialBasis,
          settings->singValThreshold, settings->lbDamping, settings->lbCycles, settings->carterCycles);
    }
  }
  else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
    auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(new FDEPotentials<SCFMode>(
        potBundle, esiPot, nonAddExchange,
        std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController())), ecpInt, pcm));
    OutputControl::nOut << "  Applying a fermi-shifted Huzinaga operator with an Fermi level of "
                        << settings->fermiShift << " Hartree." << std::endl;
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode>>(
        activeSystem, environmentSystems, *settings, superPot_AA, topDown, grid, gridCutOff, eConts, settings->fermiShift);
  }
  else {
    kin = std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController()));
  }
  if (topDown) {
    return std::shared_ptr<PotentialBundle<SCFMode>>(
        new PBEPotentials<SCFMode>(potBundle, esiPot, nonAddExchange, kin, ecpInt, pcm));
  }
  else {
    return std::shared_ptr<PotentialBundle<SCFMode>>(
        new FDEPotentials<SCFMode>(potBundle, esiPot, nonAddExchange, kin, ecpInt, pcm));
  }
  assert(false);
}

template class FDEPotentialBundleFactory<Options::SCF_MODES::RESTRICTED>;
template class FDEPotentialBundleFactory<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
