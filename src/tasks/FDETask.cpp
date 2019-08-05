/**
 * @file   FDETask.cpp
 * @author Jan Unsleber, Kevin Klahr, Thomas Dresselhaus
 *
 * @date   last reworked on Nov 11. 2016
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
#include "tasks/FDETask.h"
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "potentials/bundles/ESIPotentials.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "potentials/BUReconstructionPotential.h"
#include "potentials/HuzinagaFDEProjectionPotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/ECPInteractionPotential.h"
#include "io/FormattedOutput.h"
#include "grid/GridControllerFactory.h"
#include "data/matrices/MatrixInBasis.h"
#include "potentials/NAddFuncPotential.h"
#include "scf/Scf.h"
#include "misc/WarningTracker.h"
#include "potentials/LevelshiftHybridPotential.h"
/* Include Std and External Headers */
#include <cassert>
#include <limits>
#include <memory>


namespace Serenity {

template<Options::SCF_MODES SCFMode>
FDETask<SCFMode>::FDETask(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems):
    _activeSystem(activeSystem),
    _environmentSystems(environmentSystems),
    _supersystemgrid(nullptr),
    _finalGrid(nullptr){
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
  std::vector<std::shared_ptr<Atom> > superSystemAtoms;
  if (!_supersystemgrid) superSystemAtoms.insert(superSystemAtoms.end(),
      _activeSystem->getAtoms().begin(),
      _activeSystem->getAtoms().end());

  _activeSystem->setDiskMode(false);
  for (auto sys : _environmentSystems){
    if (sys->getSettings().scfMode == RESTRICTED) {
      sys->getElectronicStructure<RESTRICTED>();
    } else {
      sys->getElectronicStructure<UNRESTRICTED>();
    }
    sys->setDiskMode(true);
    if (!_supersystemgrid){
      double cutoff = settings.gridCutOff;
      if (cutoff<0.0) cutoff = std::numeric_limits<double>::infinity();
      for (auto atom : sys->getAtoms()){
        for (auto check : _activeSystem->getAtoms()){
          double dist = distance(*atom,*check);
          if (dist< cutoff) {
            superSystemAtoms.push_back(atom);
            break;
          }
        }
      }
    }
  }

  if (!_supersystemgrid){
    // geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtoms);
    superSystemGeometry->deleteIdenticalAtoms();
    // supersystem grid
    Options::GRID_PURPOSES gridacc = (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL:Options::GRID_PURPOSES::DEFAULT;
    _supersystemgrid = GridControllerFactory::produce(
        superSystemGeometry, _activeSystem->getSettings(), gridacc);
  }
  // Set the grid controller to the supersystem grid if "exact" methods are used.
  if(settings.embeddingMode!=Options::KIN_EMBEDDING_MODES::NADD_FUNC){
    _activeSystem->setGridController(_supersystemgrid);
  }

  // geometries of the environment subsystems
  std::vector<std::shared_ptr<const Geometry> > envGeometries;
  for (auto sys : _environmentSystems){
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
  std::vector<std::shared_ptr<EnergyComponentController> > eConts;
  double totEnvEnergy = 0.0;
  //Check if interaction energies of active systems are present at this point
  bool allInteractionsPresent=_activeSystem->template getElectronicStructure<SCFMode>()->
      checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY);
  eConts.push_back(_activeSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
  for (auto sys : _environmentSystems){
    assert((sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
        sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF)
        && "Unknown electronic structure theory!");
    auto energyContribution = (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)?
        ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY
        : ENERGY_CONTRIBUTIONS::HF_ENERGY;
    if (sys->getSettings().scfMode == SCFMode) {
      //Add isolated energies to totEnvEnergy; we need them either way.
      if (sys->template getElectronicStructure<SCFMode>()->checkEnergy(energyContribution)){
        totEnvEnergy += sys->template getElectronicStructure<SCFMode>()->
            getEnergy(energyContribution);
      }
      //Check interaction energies of every env system are present at this point
      allInteractionsPresent*=sys->template getElectronicStructure<SCFMode>()->
          checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY)
          ||(sys->template getElectronicStructure<SCFMode>()->
              checkEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY)
              &&sys->template getElectronicStructure<SCFMode>()->
              checkEnergy(ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
      eConts.push_back(sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController());
    } else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        WarningTracker::printWarning("Warning: Calculating environment system " + sys->getSettings().name +  " in restricted scfMode.\n Will proceed by assuming that alpha and beta densities are equal for this system...\n", iOOptions.printSCFCycleInfo);
        //Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->checkEnergy(energyContribution)){
          totEnvEnergy += sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
              getEnergy(energyContribution);
        }
        //Check interaction energies of every env system are present at this point
        allInteractionsPresent*=sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
            checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY)
            ||(sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
                checkEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY)
                &&sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
                checkEnergy(ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
        eConts.push_back(sys->template getElectronicStructure<RESTRICTED>()->getEnergyComponentController());
      } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        WarningTracker::printWarning("Warning: Calculating environment system " + sys->getSettings().name +  " in unrestricted scfMode.\n Will proceed by summing up alpha and beta densities for this system...\n", iOOptions.printSCFCycleInfo);
        //Add KS-DFT energies to totEnvEnergy; we need them either way.
        if (sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->checkEnergy(energyContribution)){
          totEnvEnergy += sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->
              getEnergy(energyContribution);
        }
        //Check interaction energies of every env system are present at this point
        allInteractionsPresent*=sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->
            checkEnergy(ENERGY_CONTRIBUTIONS::FDE_EMBEDDED_KS_DFT_ENERGY)
            ||(sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->
                checkEnergy(ENERGY_CONTRIBUTIONS::HF_ENERGY)
                &&sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->
                checkEnergy(ENERGY_CONTRIBUTIONS::FDE_INTERACTION_ENERGY));
        eConts.push_back(sys->template getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController());
      } else {
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
  if(allInteractionsPresent){
    for (auto sys : _environmentSystems){
      if (sys->getSettings().scfMode == SCFMode) {
        totEnvEnergy += 0.5 * sys->template getElectronicStructure<SCFMode>()->
            getEnergy(ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
      } else {
        if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->
              getEnergy(ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
        } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
          totEnvEnergy += 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->
              getEnergy(ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
        } else {
          assert(false);
        }

      }
    }
    totEnvEnergy -= 0.5 * _activeSystem->template getElectronicStructure<SCFMode>()->
        getEnergy(ENERGY_CONTRIBUTIONS::FDE_ELECTROSTATIC_INTERACTIONS);
  } else{
    WarningTracker::printWarning("WARNING: Not all interaction energies available. The environment energy is\n \
        calculated using only isolated KS-DFT energies. Be aware that this does not lead to\n \
        the correct supersystem energy. If you want the supersystem energy, please calculate\n \
        the interaction energies for every subsystem (e.g., by choosing them as active system\n \
        in a FDE or FAT calculation.)\n", iOOptions.printSCFCycleInfo);
  }


  // list of environment density matrices (their controllers)
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDensities;
  for (auto sys : _environmentSystems){
    if (sys->getSettings().scfMode == SCFMode) {
      envDensities.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    } else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        //Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        //Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin = sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        envDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      } else {
        assert(false);
      }
    }
  }

  /*
   * Create Potentials
   */
  std::shared_ptr<PotentialBundle<SCFMode> > activeSysPot;
  if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
    activeSysPot = _activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  } else if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    activeSysPot = _activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  } else {
    std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
    assert(false);
  }
  // ECP
  // Gather environment atoms
  std::vector<std::shared_ptr<Atom> > envAtoms;
  for (const auto& env : _environmentSystems) {
    for (const auto& atom : env->getGeometry()->getAtoms()) {
      envAtoms.push_back(atom);
    }
  }
  std::shared_ptr<Potential<SCFMode> > ecpInt
    (new ECPInteractionPotential<SCFMode>(
        _activeSystem,
        _activeSystem->getGeometry()->getAtoms(),
        envAtoms,
        envDensities,
        _activeSystem->getBasisController()));
  std::shared_ptr<PotentialBundle<SCFMode> > esiPot;
  if (actsettings.dft.densityFitting==Options::DENS_FITS::RI){
    std::vector<std::shared_ptr<BasisController> > envAuxBasis;
    for(auto& env : _environmentSystems){
      envAuxBasis.push_back(env->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
    }
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            _activeSystem,
            _environmentSystems,
            _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            _activeSystem->getGeometry(),
            envDensities,
            envGeometries,
            _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
            envAuxBasis));
  } else {
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            _activeSystem,
            _environmentSystems,
            _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            _activeSystem->getGeometry(),
            envDensities,
            envGeometries));
  }
  auto naddXCfunc = FunctionalClassResolver::resolveFunctional(settings.naddXCFunc);
  // Issue a warning if exact exchange should be evaluated with non-orthogonal orbitals
  // without any correction respecting this.
  if(naddXCfunc.isHybrid()
      and (settings.embeddingMode == Options::KIN_EMBEDDING_MODES::NADD_FUNC
          ||settings.embeddingMode == Options::KIN_EMBEDDING_MODES::NONE)) {
    WarningTracker::printWarning("Warning: Exact exchange is calculated with orbitals that may not be orthogonal.",iOOptions.printSCFCycleInfo);
  }

  auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode> >(
      _activeSystem,
      _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      envDensities,
      _supersystemgrid,
      naddXCfunc,
      std::make_pair(true,(settings.gridCutOff<0)?eConts:std::vector<std::shared_ptr<EnergyComponentController> >(0)),
      (settings.gridCutOff<0.0)?true:false);

  std::shared_ptr<Potential<SCFMode> > kin;
  if(settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NADD_FUNC){
    kin = std::make_shared<NAddFuncPotential<SCFMode> >(
        _activeSystem,
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        _supersystemgrid,
        FunctionalClassResolver::resolveFunctional(settings.naddKinFunc),
        std::make_pair(false,(settings.gridCutOff<0)?eConts:std::vector<std::shared_ptr<EnergyComponentController> >(0)),
        (settings.gridCutOff<0.0)?true:false);
  }else if(settings.embeddingMode==Options::KIN_EMBEDDING_MODES::LEVELSHIFT){
    kin = std::make_shared<LevelshiftHybridPotential<SCFMode> >(
        _activeSystem,
        _environmentSystems,
        settings.levelShiftParameter,
        (settings.basisFunctionRatio == 0.0)? false : true,
         settings.basisFunctionRatio,
         settings.borderAtomThreshold,
         _supersystemgrid,
         FunctionalClassResolver::resolveFunctional(settings.naddKinFunc),
         false,
         std::make_pair(false,(settings.gridCutOff<0)?eConts:std::vector<std::shared_ptr<EnergyComponentController> >(0)));
  }else if(settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HUZINAGA){
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode> >(
        _activeSystem,
        _environmentSystems,
        naddXCfunc,
        settings.truncateProjector,
        settings.projecTruncThresh,
        false,
        false,
        settings.distantKinFunc,
        settings.basisFunctionRatio,
        settings.borderAtomThreshold,
        _supersystemgrid,
        settings.naddKinFunc,
        settings.gridCutOff,
        eConts);
  }else if(settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HOFFMANN){
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode> >(
        _activeSystem,
        _environmentSystems,
        naddXCfunc,
        settings.truncateProjector,
        settings.projecTruncThresh,
        true,
        false,
        settings.distantKinFunc,
        settings.basisFunctionRatio,
        settings.borderAtomThreshold,
        _supersystemgrid,
        settings.naddKinFunc,
        settings.gridCutOff,
        eConts);
  }else if(settings.embeddingMode==Options::KIN_EMBEDDING_MODES::RECONSTRUCTION){
    kin=std::make_shared<BUReconstructionPotential<SCFMode>>(
        _activeSystem,
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        _supersystemgrid,
        naddXCfunc,
        _environmentSystems,
        settings.smoothFactor,
        settings.potentialBasis.empty()? _activeSystem->getBasisController()->getBasisString() : settings.potentialBasis,
            settings.singValThreshold,
            settings.lbDamping,
            settings.lbCycles,
            settings.carterCycles);
  }else{
    kin = std::shared_ptr<Potential<SCFMode> >(
        new ZeroPotential<SCFMode>(_activeSystem->getBasisController()));
  }

  /*
   * Bundle potentials
   */
  auto fdePot = std::shared_ptr<PotentialBundle<SCFMode> >(
      new FDEPotentials<SCFMode>(
          activeSysPot,
          esiPot,
          naddXCPot,
          kin,
          ecpInt));

  /*
   * non-additive dispersion correction
   */
  auto nadDispCorrection = 0.0;
  auto supersystemGeometry = std::make_shared<Geometry>();
  *supersystemGeometry += *_activeSystem->getGeometry();
  auto actDispCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
      _activeSystem->getGeometry(),settings.naddXCFunc);
  nadDispCorrection -= actDispCorrection;
  for (auto sys : _environmentSystems){
    *supersystemGeometry += *sys->getGeometry();
    nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
        sys->getGeometry(),settings.naddXCFunc);
  }
  nadDispCorrection += DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
      supersystemGeometry,settings.naddXCFunc);


  /*
   * Attach everything and run scf
   */

  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION,actDispCorrection);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP,nadDispCorrection);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES,totEnvEnergy);
  Scf<SCFMode>::perform(actsettings,es,fdePot);
  if(settings.gridCutOff>0.0 and settings.finalGrid and settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NADD_FUNC){
    if (!_finalGrid){
      std::cout << "Creating final grid for energy evaluation" << std::endl;
      // atoms of all subsystems
      auto finalGridGeometry=std::make_shared<Geometry>();
      *finalGridGeometry+=*_activeSystem->getGeometry();

      for (auto sys : _environmentSystems){
        *finalGridGeometry+=*sys->getGeometry();
      }

      // supersystem grid
      Options::GRID_PURPOSES finalGridacc = (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL:Options::GRID_PURPOSES::DEFAULT;
      _finalGrid = GridControllerFactory::produce(
          finalGridGeometry, _activeSystem->getSettings(), finalGridacc);
    }
    naddXCPot = std::make_shared<NAddFuncPotential<SCFMode> >(
        _activeSystem,
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        _finalGrid,
        naddXCfunc,
        std::make_pair(true,eConts),
        true);
    kin = std::make_shared<NAddFuncPotential<SCFMode> >(
        _activeSystem,
        _activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        _finalGrid,
        FunctionalClassResolver::resolveFunctional(settings.naddKinFunc),
        std::make_pair(false,eConts),
        true);

    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC,naddXCPot->getEnergy(_activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix()));
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC,kin->getEnergy(_activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix()));
  }
  /*
   * Print final results if sensible
   */
  if (iOOptions.printSCFResults and not (settings.gridCutOff>0.0 and !settings.finalGrid)){
    printSubSectionTitle("Final SCF Results");
    eCont->printAllComponents();
    if (iOOptions.printFinalOrbitalEnergies){
      print("");
      printSmallCaption("Orbital Energies");
      es->printMOEnergies();
      print("");
    }
  }

}

template class FDETask<Options::SCF_MODES::RESTRICTED>;
template class FDETask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
