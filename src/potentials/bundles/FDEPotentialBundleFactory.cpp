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
#include "basis/BasisFunctionMapper.h"
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
#include "potentials/HuzinagaProjectionPotential.h"
#include "potentials/LevelshiftHybridPotential.h"
#include "potentials/LoewdinFDEProjectionPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/PCMPotential.h"
#include "potentials/TDReconstructionPotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/ALMOPotentials.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/bundles/PBEPotentials.h"
namespace Serenity {

template<>
std::shared_ptr<PotentialBundle<RESTRICTED>> FDEPotentialBundleFactory<RESTRICTED>::produce(
    std::shared_ptr<SystemController> activeSystem, std::shared_ptr<DensityMatrixController<RESTRICTED>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<RESTRICTED>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
    std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex) {
  if (!_restrictedFDEFactory)
    _restrictedFDEFactory.reset(new FDEPotentialBundleFactory<RESTRICTED>);
  if (environmentSystems.size() == 0) {
    if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      return activeSystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
    }
    else if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      return activeSystem->getPotentials<RESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    }
  }
  return _restrictedFDEFactory->getOrProduce(activeSystem, activeDensMatController, environmentSystems,
                                             envDensMatController, settings, grid, supersystem, topDown, noSuperRecon,
                                             gridCutOff, eConts, firstPassiveSystemIndex);
}

template<>
std::shared_ptr<PotentialBundle<UNRESTRICTED>> FDEPotentialBundleFactory<UNRESTRICTED>::produce(
    std::shared_ptr<SystemController> activeSystem,
    std::shared_ptr<DensityMatrixController<UNRESTRICTED>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<UNRESTRICTED>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
    std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex) {
  if (!_unrestrictedFDEFactory)
    _unrestrictedFDEFactory.reset(new FDEPotentialBundleFactory<UNRESTRICTED>);
  if (environmentSystems.size() == 0) {
    if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      return activeSystem->getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
    }
    else if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      return activeSystem->getPotentials<UNRESTRICTED, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    }
  }
  return _unrestrictedFDEFactory->getOrProduce(activeSystem, activeDensMatController, environmentSystems,
                                               envDensMatController, settings, grid, supersystem, topDown, noSuperRecon,
                                               gridCutOff, eConts, firstPassiveSystemIndex);
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<PotentialBundle<SCFMode>> FDEPotentialBundleFactory<SCFMode>::produceNew(
    std::shared_ptr<SystemController> activeSystem, std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
    std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex) {
  if (environmentSystems.size() != envDensMatController.size())
    throw SerenityError("ERROR: Inconsistent number of environment density matrices and environment systems!");

  std::shared_ptr<PotentialBundle<SCFMode>> potBundle;
  if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    potBundle = activeSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  }
  else if (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    potBundle = activeSystem->getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  }
  else {
    throw SerenityError("Nonexistent electronicStructureTheory requested. Options are HF and DFT.");
  }
  auto densFitJ = activeSystem->getSettings().basis.densFitJ;
  // The collected atoms and geometries will be needed further down for ECPs and Coulomb interaction.
  std::vector<std::shared_ptr<const Geometry>> environmentGeometries;
  std::vector<std::shared_ptr<BasisController>> environmentAuxiliaryBasisSets;
  std::vector<std::shared_ptr<Atom>> environmentAtoms;
  for (auto env : environmentSystems) {
    environmentGeometries.push_back(env->getGeometry());
  }
  for (auto envGeom : environmentGeometries)
    environmentAtoms.insert(environmentAtoms.begin(), envGeom->getAtoms().begin(), envGeom->getAtoms().end());
  std::shared_ptr<PotentialBundle<SCFMode>> esiPot;
  if (densFitJ != Options::DENS_FITS::NONE) {
    for (auto env : environmentSystems) {
      environmentAuxiliaryBasisSets.push_back(env->getAuxBasisController(Options::AUX_BASIS_PURPOSES::COULOMB, densFitJ));
    }
    esiPot = std::make_shared<ESIPotentials<SCFMode>>(
        activeSystem, environmentSystems, activeDensMatController, activeSystem->getGeometry(), envDensMatController,
        environmentGeometries, activeSystem->getAuxBasisController(Options::AUX_BASIS_PURPOSES::COULOMB, densFitJ),
        environmentAuxiliaryBasisSets, firstPassiveSystemIndex, settings->partialChargesForCoulombInt, settings->chargeModel);
  }
  else {
    esiPot = std::make_shared<ESIPotentials<SCFMode>>(
        activeSystem, environmentSystems, activeDensMatController, activeSystem->getGeometry(), envDensMatController,
        environmentGeometries, firstPassiveSystemIndex, settings->partialChargesForCoulombInt, settings->chargeModel);
  }
  auto ecpInt_total =
      std::make_shared<ECPInteractionPotential<SCFMode>>(activeSystem, activeSystem->getGeometry()->getAtoms(), environmentAtoms,
                                                         envDensMatController, activeSystem->getBasisController());

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
  // TODO add a decomposition of the cavity formation energy which is suitable for FDE.
  auto pcm = std::make_shared<PCMPotential<SCFMode>>(
      settings->pcm, activeSystem->getBasisController(), activeSystem->getGeometry(),
      (usesPCM) ? activeSystem->getMolecularSurface(MOLECULAR_SURFACE_TYPES::FDE) : nullptr, nullptr,
      (usesPCM) ? activeSystem->getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::FDE)
                : nullptr,
      envElecPots);

  if (settings->naddXCFuncList.size() > 0) {
    return buildMixedEmbedding(activeSystem, activeDensMatController, environmentSystems, envDensMatController,
                               settings, grid, topDown, gridCutOff, eConts, potBundle, esiPot, pcm, ecpInt_total);
  }
  else {
    auto naddXCfunc = settings->customNaddXCFunc.basicFunctionals.size() ? Functional(settings->customNaddXCFunc)
                                                                         : resolveFunctional(settings->naddXCFunc);
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
        (topDown || gridCutOff < 0.0),
        (settings->embeddingMode != Options::KIN_EMBEDDING_MODES::RECONSTRUCTION && gridCutOff < 0.0));

    std::shared_ptr<Potential<SCFMode>> kin;
    if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
      kin = std::make_shared<NAddFuncPotential<SCFMode>>(
          activeSystem, activeDensMatController, envDensMatController, grid,
          settings->customNaddKinFunc.basicFunctionals.size() ? Functional(settings->customNaddKinFunc)
                                                              : resolveFunctional(settings->naddKinFunc),
          std::make_pair(false, (gridCutOff < 0) ? eConts : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
          (gridCutOff < 0.0));
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
      kin = std::make_shared<LevelshiftHybridPotential<SCFMode>>(
          activeSystem, environmentSystems, settings->levelShiftParameter, (settings->basisFunctionRatio == 0.0) ? false : true,
          settings->basisFunctionRatio, settings->borderAtomThreshold, grid,
          settings->customLongRangeNaddKinFunc.basicFunctionals.size() ? Functional(settings->customLongRangeNaddKinFunc)
                                                                       : resolveFunctional(settings->longRangeNaddKinFunc),
          topDown,
          std::make_pair(false, (gridCutOff < 0 && !topDown) ? eConts
                                                             : std::vector<std::shared_ptr<EnergyComponentController>>(0)));
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
      auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(new FDEPotentials<SCFMode>(
          potBundle, esiPot, {nonAddExchange},
          {std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController()))},
          ecpInt_total, pcm));
      kin = std::make_shared<HuzinagaProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings,
                                                                   superPot_AA, topDown, grid, gridCutOff, eConts);
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::LOEWDIN) {
      if (topDown)
        throw SerenityError("The Loewdin algorithm is not intended to be used in the top-down scenario!");
      kin = std::make_shared<LoewdinFDEProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings);
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
      auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(new FDEPotentials<SCFMode>(
          potBundle, esiPot, {nonAddExchange},
          {std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController()))},
          ecpInt_total, pcm));
      kin = std::make_shared<HoffmannProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings,
                                                                   superPot_AA, topDown, grid, gridCutOff, eConts);
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::ALMO) {
      for (auto envSys : environmentSystems) {
        if (activeSystem->getSettings().dft.functional != envSys->getSettings().dft.functional)
          throw SerenityError("The same exchange-correlation functional is required for all subsystems using ALMOs.");
      }
      if (settings->naddXCFunc != activeSystem->getSettings().dft.functional)
        WarningTracker::printWarning("Warning: The non-additive exchange-correlation functional is ignored for ALMOs.",
                                     iOOptions.printSCFCycleInfo);
      return std::unique_ptr<PotentialBundle<SCFMode>>(new ALMOPotentials<SCFMode>(activeSystem, environmentSystems));
    }
    else if (settings->embeddingMode == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION) {
      if (topDown) {
        kin = std::make_shared<TDReconstructionPotential<SCFMode>>(
            activeSystem, supersystem, environmentSystems, settings->smoothFactor,
            settings->potentialBasis.empty() ? environmentSystems[0]->getBasisController()->getBasisString()
                                             : settings->potentialBasis,
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
          potBundle, esiPot, {nonAddExchange},
          {std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController()))},
          ecpInt_total, pcm));
      OutputControl::nOut << "  Applying a fermi-shifted Huzinaga operator with a fermi level of "
                          << settings->fermiShift << " Hartree." << std::endl;
      kin = std::make_shared<HuzinagaProjectionPotential<SCFMode>>(activeSystem, environmentSystems, *settings, superPot_AA,
                                                                   topDown, grid, gridCutOff, eConts, settings->fermiShift);
    }
    else {
      kin = std::shared_ptr<Potential<SCFMode>>(new ZeroPotential<SCFMode>(activeSystem->getBasisController()));
    }
    return std::unique_ptr<PotentialBundle<SCFMode>>(
        new FDEPotentials<SCFMode>(potBundle, esiPot, {nonAddExchange}, {kin}, ecpInt_total, pcm));
    assert(false);
  }
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<PotentialBundle<SCFMode>> FDEPotentialBundleFactory<SCFMode>::buildMixedEmbedding(
    std::shared_ptr<SystemController> activeSystem, std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
    std::vector<std::shared_ptr<SystemController>> environmentSystems,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
    const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid, bool topDown,
    double gridCutOff, std::vector<std::shared_ptr<EnergyComponentController>> eConts,
    std::shared_ptr<PotentialBundle<SCFMode>> potBundle, std::shared_ptr<PotentialBundle<SCFMode>> esiPot,
    std::shared_ptr<PCMPotential<SCFMode>> pcm, std::shared_ptr<ECPInteractionPotential<SCFMode>> ecpInt_total) {
  if (environmentSystems.size() + 1 != settings->embeddingModeList.size()) {
    throw SerenityError("You want to use list input for embedding; You need to specify the same number of embedding "
                        "modes as environment and active subsystems!");
  }
  for (auto embMode : settings->embeddingModeList) {
    if (embMode == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION)
      throw SerenityError("Potential Reconstruction in combined Approximate/Exact embedding not supported");
    if (embMode == Options::KIN_EMBEDDING_MODES::ALMO)
      throw SerenityError("ALMOs in combined Approximate/Exact embedding not supported");
  }
  if (envDensMatController.size() != environmentSystems.size()) {
    throw SerenityError("Different number of environment systems and density matrices!");
  }
  if (settings->naddXCFuncList.size() != 2) {
    throw SerenityError("You need to specify 2 nadd-xc-functionals for combined approx/exact calculation!");
  }
  // Variables for list input
  auto naddXCFuncs = settings->naddXCFuncList;
  auto embeddingModes = settings->embeddingModeList;
  // Variables
  std::vector<std::shared_ptr<NAddFuncPotential<SCFMode>>> nonAddXCVector;
  std::vector<std::shared_ptr<Potential<SCFMode>>> nonAddKinVector;
  BasisFunctionMapper basisFunctionMapper(activeSystem->getBasisController());
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> _BtoAProjections;
  std::vector<std::shared_ptr<SystemController>> exactlyTreatedSystems;
  std::vector<std::shared_ptr<SystemController>> otherSystems;
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> exactlyTreatedenvDensMatController;
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> otherDensMatController;

  std::vector<std::shared_ptr<EnergyComponentController>> eConts_Sub_Exact;
  std::vector<std::shared_ptr<EnergyComponentController>> eConts_Sub_Approx;

  if (eConts.size() != 0) {
    std::vector<std::shared_ptr<EnergyComponentController>> eConts_Sub_Exact = {eConts[0]};
    std::vector<std::shared_ptr<EnergyComponentController>> eConts_Sub_Approx = {eConts[0]};
  }

  // Seperate between exactly and approximately treated environment systems
  for (unsigned int iEnv = 0; iEnv < environmentSystems.size(); ++iEnv) {
    if (embeddingModes[iEnv + 1] == Options::KIN_EMBEDDING_MODES::HUZINAGA ||
        embeddingModes[iEnv + 1] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
      exactlyTreatedSystems.push_back(environmentSystems[iEnv]);
      exactlyTreatedenvDensMatController.push_back(envDensMatController[iEnv]);
      if (!eConts_Sub_Exact.empty())
        eConts_Sub_Exact.push_back(eConts[iEnv + 1]);
      auto basisControllerB = environmentSystems[iEnv]->getBasisController();
      _BtoAProjections.push_back(basisFunctionMapper.getSparseProjection(basisControllerB));
    }
    else {
      otherSystems.push_back(environmentSystems[iEnv]);
      otherDensMatController.push_back(envDensMatController[iEnv]);
      if (!eConts_Sub_Approx.empty())
        eConts_Sub_Approx.push_back(eConts[iEnv + 1]);
    }
  }
  // Some debug output
  OutputControl::dOut << "Exactly treated subsystems - " << exactlyTreatedSystems.size() + 1 << std::endl;
  OutputControl::dOut << "Approximately treated subsystems - " << otherSystems.size() << std::endl;

  // NaddXCPotential for exactly treated systems
  auto naddXCfunc = resolveFunctional(naddXCFuncs[0]);
  nonAddXCVector.push_back(std::make_shared<NAddFuncPotential<SCFMode>>(
      activeSystem, activeDensMatController, exactlyTreatedenvDensMatController, grid, naddXCfunc,
      std::make_pair(true, (gridCutOff < 0) ? eConts_Sub_Exact : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
      (topDown || gridCutOff < 0.0),
      (settings->embeddingMode != Options::KIN_EMBEDDING_MODES::RECONSTRUCTION && gridCutOff < 0.0)));

  // NaddXCPotential for not exactly treated subsystems
  naddXCfunc = resolveFunctional(naddXCFuncs[1]);
  nonAddXCVector.push_back(std::make_shared<NAddFuncPotential<SCFMode>>(
      activeSystem, activeDensMatController, exactlyTreatedenvDensMatController, _BtoAProjections,
      otherDensMatController, grid, naddXCfunc,
      std::make_pair(true, (gridCutOff < 0) ? eConts_Sub_Approx : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
      (topDown || gridCutOff < 0.0),
      (settings->embeddingMode != Options::KIN_EMBEDDING_MODES::RECONSTRUCTION && gridCutOff < 0.0)));
  // NaddKin
  auto naddKinfunc = settings->customNaddKinFunc.basicFunctionals.size() ? Functional(settings->customNaddKinFunc)
                                                                         : resolveFunctional(settings->naddKinFunc);
  nonAddKinVector.push_back(std::make_shared<NAddFuncPotential<SCFMode>>(
      activeSystem, activeDensMatController, exactlyTreatedenvDensMatController, _BtoAProjections,
      otherDensMatController, grid, naddKinfunc,
      std::make_pair(false, (gridCutOff < 0) ? eConts_Sub_Approx : std::vector<std::shared_ptr<EnergyComponentController>>(0)),
      (gridCutOff < 0.0)));

  // Projector for exactly treated systems
  if (embeddingModes[0] == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
    // Set up environment geometry and list of environment geometries for different potentials
    std::shared_ptr<Geometry> envGeomExactCombined = std::make_shared<Geometry>();
    std::vector<std::shared_ptr<const Geometry>> envGeometry;
    for (auto sys : exactlyTreatedSystems) {
      (*envGeomExactCombined) += (*sys->getGeometry());
      envGeometry.push_back(sys->getGeometry());
    }
    envGeomExactCombined->deleteIdenticalAtoms();
    auto superPot_AA = std::shared_ptr<PotentialBundle<SCFMode>>(
        new FDEPotentials<SCFMode>(potBundle, esiPot, nonAddXCVector, nonAddKinVector, ecpInt_total, pcm));
    settings->naddXCFunc = naddXCFuncs[0];
    nonAddKinVector.push_back(std::make_shared<HuzinagaProjectionPotential<SCFMode>>(
        activeSystem, exactlyTreatedSystems, *settings, superPot_AA, topDown, grid, gridCutOff, eConts_Sub_Exact));
  }
  else if (embeddingModes[0] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
    nonAddKinVector.push_back(std::make_shared<LevelshiftHybridPotential<SCFMode>>(
        activeSystem, exactlyTreatedSystems, settings->levelShiftParameter,
        (settings->basisFunctionRatio == 0.0) ? false : true, settings->basisFunctionRatio, settings->borderAtomThreshold, grid,
        settings->customLongRangeNaddKinFunc.basicFunctionals.size() ? Functional(settings->customLongRangeNaddKinFunc)
                                                                     : resolveFunctional(settings->longRangeNaddKinFunc),
        topDown,
        std::make_pair(false, (gridCutOff < 0 && !topDown) ? eConts_Sub_Exact
                                                           : std::vector<std::shared_ptr<EnergyComponentController>>(0))));
  }
  else {
    throw SerenityError("Projector not supported for combined (exact/approx.) embedding strategy!");
  }

  OutputControl::dOut << "Number of non-additive kinetic interactions - " << nonAddKinVector.size() << std::endl;
  OutputControl::dOut << "Number of non-additive exchange-correlation interactions - " << nonAddXCVector.size() << std::endl;

  return std::unique_ptr<PotentialBundle<SCFMode>>(
      new FDEPotentials<SCFMode>(potBundle, esiPot, nonAddXCVector, nonAddKinVector, ecpInt_total, pcm));
}

template class FDEPotentialBundleFactory<Options::SCF_MODES::RESTRICTED>;
template class FDEPotentialBundleFactory<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
