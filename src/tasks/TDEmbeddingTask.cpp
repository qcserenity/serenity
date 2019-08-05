/**
 * @file   TDEmbeddingTask.cpp
 *
 * @date   Apr 23, 2014
 * @author Jan Unsleber, Thomas Dresselhaus
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
#include "tasks/TDEmbeddingTask.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/OrbitalController.h"
#include "data/matrices/FockMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "grid/GridControllerFactory.h"
#include "integrals/OneElectronIntegralController.h"
#include "misc/SystemSplittingTools.h"
#include "misc/WarningTracker.h"
#include "postHF/MPn/MP2.h"
#include "postHF/MPn/RIMP2.h"
#include "scf/Scf.h"
#include "system/SystemController.h"
/* Tasks */
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
/* Potentials and bundles*/
#include "potentials/bundles/PBEPotentials.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/ECPInteractionPotential.h"
#include "potentials/HuzinagaFDEProjectionPotential.h"
#include "potentials/LevelshiftHybridPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/TDReconstructionPotential.h"
#include "potentials/ZeroPotential.h"

/* Include Std and External Headers */
#include <memory>

namespace Serenity {
using namespace std;
template<Options::SCF_MODES SCFMode>
TDEmbeddingTask<SCFMode>::TDEmbeddingTask(
    std::shared_ptr<SystemController> activeSystem,
    std::shared_ptr<SystemController> environmentSystem):
    _activeSystem(activeSystem),
    _environmentSystem(environmentSystem){
}

template<Options::SCF_MODES SCFMode>
void TDEmbeddingTask<SCFMode>::run() {
  // Check input for obvious misunderstanding of the keywords.
  checkInput();
  /*
   * Trigger basis generation to add basis functions to
   * active and environment atoms
   */
  if(!settings.useEnvSys){
    _activeSystem->setBasisController(nullptr);
    _environmentSystem->setBasisController(nullptr);
    _activeSystem->getAtomCenteredBasisController()->getBasis();
    _environmentSystem->getAtomCenteredBasisController()->getBasis();
  }
  /*
   * Load or create a supersystem from active and environment systems.
   * This may trigger a supersystem SCF.
   */
  std::shared_ptr<SystemController> supersystem = setUpSupersystem();
  /*
   * Funny enough the "supersystemGrid" is build with the settings of the active system while the
   * grid of the supersystem is build with the settings of the environment.
   */
  std::shared_ptr<AtomCenteredGridController> supersystemGrid=_activeSystem->getAtomCenteredGridController();
  /*
   * Partitioning
   * Supersystem orbital localization to atoms
   */
  LocalizationTask superSystemOrbLocalization(supersystem);
  superSystemOrbLocalization.settings.locType = settings.locType;
  superSystemOrbLocalization.run();

  auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
      supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients(),
      supersystem->getOneElectronIntegralController()->getOverlapIntegrals(),
      supersystem->getAtomCenteredBasisController()->getBasisIndices());

  // Print some information about the picked orbitals and
  //   the level shift parameter
  if (settings.embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT){
    std::cout <<std::endl;
    std::cout<< "Levelshifting/projectiong active orbitals by: "<<settings.levelShiftParameter<<" Hartree."<<std::endl;
  }
  std::cout <<std::endl;
  auto listOfActiveOrbs = SystemSplittingTools<SCFMode>::partitionOrbitals(
      orbitalPopulations,
      supersystem,
      _activeAtoms,
      settings.orbitalThreshold,
      _activeSystem->getNOccupiedOrbitals<SCFMode>(),
      settings.enforceCharges);

  const auto suffix=getDensityNameSuffix();
  SpinPolarizedData<SCFMode, unsigned int> nActiveOrbs(0);
  for_spin(listOfActiveOrbs,suffix,nActiveOrbs){
    std::cout << std::endl << "Choosing " << suffix_spin << " orbitals for active system"  << std::endl;
    for (unsigned int i=0; i<listOfActiveOrbs_spin.size();++i){
      if (listOfActiveOrbs_spin[i]){
        std::cout << i+1 << " ";
        nActiveOrbs_spin += 1;
      }
      if ((i+1)%10 == 0) std::cout << std::endl;
    }
    std::cout << std::endl << "Number of active " << suffix_spin <<  " orbitals: "<< nActiveOrbs_spin << std::endl;
    std::cout << std::endl;
  };


  /* Split the electronic structure of the supersystem */
  auto electronicStructurePair = SystemSplittingTools<SCFMode>::splitElectronicStructure(
      supersystem,listOfActiveOrbs);
  auto& actES = electronicStructurePair.first;
  auto& envES = electronicStructurePair.second;
  if (settings.useEnvSys) envES = _environmentSystem->template getElectronicStructure<SCFMode>();
  _activeSystem->getGeometry()->addAsDummy(*(_environmentSystem->getGeometry()));
  _environmentSystem->getGeometry()->addAsDummy(*(_activeSystem->getGeometry()));
  _activeSystem->getGeometry()->deleteIdenticalAtoms();
  _environmentSystem->getGeometry()->deleteIdenticalAtoms();
  _environmentSystem->setBasisController(dynamic_pointer_cast<AtomCenteredBasisController>(
      supersystem->getAtomCenteredBasisController()));
  _activeSystem->setBasisController(dynamic_pointer_cast<AtomCenteredBasisController>(
      supersystem->getAtomCenteredBasisController()));
  _activeSystem->setElectronicStructure<SCFMode>(actES);
  _environmentSystem->setElectronicStructure<SCFMode>(envES);
  _environmentSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->updateDensityMatrix();

  /*
   * Adjust the charges of the given systems
   *
   * Even though the user specifies charges the actual charges will be determined
   * by counting the occupied orbitals localized to the fragment.
   */
  SpinPolarizedData<SCFMode,unsigned int> activeSystemElectrons = SystemSplittingTools<SCFMode>::getNElectrons(
      listOfActiveOrbs).first;
  if (!settings.enforceCharges){
    // Adjust the number of electron is both systems according to
    //   the number of orbitals that have been picked, because
    //   they were localized on the active system.
    unsigned int tmpSys = 0;
    unsigned int tmpAct = 0;
    auto oldNElActive  = _activeSystem->getNElectrons<SCFMode>();
    for_spin(oldNElActive,activeSystemElectrons){
      tmpSys += oldNElActive_spin;
      tmpAct += activeSystemElectrons_spin;
    };
    int actSpin = SystemSplittingTools<SCFMode>::getSpin(activeSystemElectrons);
    _activeSystem->setSpin(actSpin);
    _activeSystem->setCharge(_activeSystem->getCharge()+ tmpSys - tmpAct);
    _environmentSystem->setSpin(supersystem->getSpin()-actSpin);
    _environmentSystem->setCharge(_environmentSystem->getCharge()- tmpSys + tmpAct);
    std::cout << "Active System Charge: " << _activeSystem->getCharge() << std::endl;
    std::cout << "Active System Spin:   " << _activeSystem->getSpin() << std::endl;
    std::cout << "Environment Charge:   " << _environmentSystem->getCharge() << std::endl;
    std::cout << "Environment Spin:     " << _environmentSystem->getSpin() << std::endl;
    std::cout <<  std::endl;
  }
  //Write supersystem to disk.
  supersystem->setDiskMode(true);
  auto densityMatrixActSys = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix();
  auto densityMatrixEnvironment =_environmentSystem->getElectronicStructure<SCFMode>()->getDensityMatrix();

  /* ========================================================= */
  /*                    Basis set truncation                   */
  /* ========================================================= */

  /*
   * The supersystem basis in which densityMatrixActSys is expressed was generated using the label
   * of the environment system. The following lines of code will:
   * 1. Set up a supersystem basis using the label of the active system
   * 2. Project densityMatrixActSys into that basis
   * 3. Truncate this basis to the most important shells.
   */

  /*
   * Set up new basis using the active systems label
   */
  auto actSysBasisController=AtomCenteredBasisControllerFactory::produce(
      _activeSystem->getGeometry(),
      _activeSystem->getSettings().basis.basisLibPath,
      _activeSystem->getSettings().basis.makeSphericalBasis,
      false,
      _activeSystem->getSettings().basis.firstECP,
      _activeSystem->getSettings().basis.label);

  _activeSystem->setBasisController(actSysBasisController);

  /*
   * Project the active systems density into the new basis
   */
  DensityMatrix<SCFMode> densInActSysBasis=SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(
        densityMatrixActSys,
        _activeSystem->getBasisController(),
        std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > (
            _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()));
  actES = std::make_shared<ElectronicStructure<SCFMode> >(
      _activeSystem->getBasisController(),
      _activeSystem->getGeometry(),
      _activeSystem->getNOccupiedOrbitals<SCFMode>());
  actES->getDensityMatrixController()->setDensityMatrix(densInActSysBasis);
  _activeSystem->setElectronicStructure<SCFMode>(actES);

  /*
   * Do the truncation
   */
  this->setActiveSystemBasis();
  _activeSystem->setGridController(supersystemGrid);
  _activeSystem->getAtomCenteredBasisController()->toHDF5(_activeSystem->getHDF5BaseName(),_activeSystem->getSettings().identifier);

  /*
   * Project active systems density matrix into new basis and attach
   */
  DensityMatrix<SCFMode> projectedDens=SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(
      densityMatrixActSys,
      _activeSystem->getBasisController(),
      std::make_shared<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > (
          _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals()));
  actES = std::make_shared<ElectronicStructure<SCFMode> >(
      _activeSystem->getBasisController(),
      _activeSystem->getGeometry(),
      _activeSystem->getNOccupiedOrbitals<SCFMode>());
  actES->getDensityMatrixController()->setDensityMatrix(projectedDens);
  _activeSystem->setElectronicStructure<SCFMode>(actES);

  // Building the density matrix controllers needed for the potential evaluation
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDensMats = {
      std::make_shared<DensityMatrixController<SCFMode> >(densityMatrixEnvironment)};
  /*
   * Create Potentials
   */
  const auto& actsettings = _activeSystem->getSettings();
  std::shared_ptr<PotentialBundle<SCFMode> > activeSysPot;
  if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
    activeSysPot = _activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  } else if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    activeSysPot = _activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  } else {
    std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
    assert(false);
  }
  std::vector<std::shared_ptr<const Geometry> > envGeometries = {_environmentSystem->getGeometry()};
  std::shared_ptr<PotentialBundle<SCFMode> > esiPot;
  if (actsettings.dft.densityFitting!=Options::DENS_FITS::RI){
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            _activeSystem,
            {_environmentSystem},
            _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            _activeSystem->getGeometry(),
            envDensMats,
            envGeometries));
  } else {
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            _activeSystem,
            {_environmentSystem},
            _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            _activeSystem->getGeometry(),
            envDensMats,
            envGeometries,
            _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
            {_environmentSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)}));
  }

  auto naddXCPot = std::shared_ptr<NAddFuncPotential<SCFMode> >( new NAddFuncPotential<SCFMode>(
      _activeSystem,
      _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      envDensMats,
      supersystemGrid,
      FunctionalClassResolver::resolveFunctional(settings.naddXCFunc),
      {true , std::vector<std::shared_ptr<EnergyComponentController> >(0)},
      true,
      settings.embeddingMode!=Options::KIN_EMBEDDING_MODES::RECONSTRUCTION));

  std::shared_ptr<Potential<SCFMode> > kin;
  if (settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NONE){
    kin = std::shared_ptr<Potential<SCFMode> >(
        new ZeroPotential<SCFMode>(_activeSystem->getBasisController()));
  } else if (settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
    kin = std::shared_ptr<Potential<SCFMode> >(
        new NAddFuncPotential<SCFMode>(
            _activeSystem,
            _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            envDensMats,
            supersystemGrid,
            FunctionalClassResolver::resolveFunctional(settings.naddKinFunc)));
  } else if (settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HUZINAGA) {
    kin = std::shared_ptr<Potential<SCFMode> >(
        new HuzinagaFDEProjectionPotential<SCFMode> (
          _activeSystem,
          {_environmentSystem},
          FunctionalClassResolver::resolveFunctional(settings.naddXCFunc),
          false,
          0.0,
          false,
          true,
          settings.nonOrthogonalCrit != Options::NON_ORTHOGONAL_CRITERION::NONE
            && settings.basisFunctionRatio != 0.0,
          settings.basisFunctionRatio,
          settings.borderAtomThreshold,
          supersystemGrid,
          settings.longRangeNaddKinFunc));
  } else if (settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HOFFMANN) {
    kin = std::shared_ptr<Potential<SCFMode> >(
        new HuzinagaFDEProjectionPotential<SCFMode> (
          _activeSystem,
          {_environmentSystem},
          FunctionalClassResolver::resolveFunctional(settings.naddXCFunc),
          false,
          0.0,
          true,
          true,
          settings.nonOrthogonalCrit != Options::NON_ORTHOGONAL_CRITERION::NONE
            && settings.basisFunctionRatio != 0.0,
          settings.basisFunctionRatio,
          settings.borderAtomThreshold,
          supersystemGrid,
          settings.longRangeNaddKinFunc));
  }else if (settings.embeddingMode==Options::KIN_EMBEDDING_MODES::RECONSTRUCTION){
    kin = std::shared_ptr<Potential<SCFMode> >(
        new TDReconstructionPotential<SCFMode>(
            _activeSystem,
            supersystem,
            {_environmentSystem},
            settings.smoothFactor,
            settings.potentialBasis.empty()?
                _environmentSystem->getBasisController()->getBasisString() : settings.potentialBasis,
                settings.singValThreshold,
                settings.lbDamping,
                settings.lbCycles,
                settings.carterCycles,
                settings.noSupRec));
  } else {
    kin = std::shared_ptr<Potential<SCFMode> >(
        new LevelshiftHybridPotential<SCFMode>(
            _activeSystem,
            {_environmentSystem},
            settings.levelShiftParameter,
            settings.nonOrthogonalCrit != Options::NON_ORTHOGONAL_CRITERION::NONE,
            settings.basisFunctionRatio,
            settings.borderAtomThreshold,
            supersystemGrid,
            FunctionalClassResolver::resolveFunctional(settings.longRangeNaddKinFunc),
            true));
  }

  // ECP TODO: Check consistency!
  std::shared_ptr<Potential<SCFMode> > ecpInt
  (new ECPInteractionPotential<SCFMode>(
      _activeSystem,
      _activeSystem->getGeometry()->getAtoms(),
      _environmentSystem->getGeometry()->getAtoms(),
      envDensMats,
      _activeSystem->getBasisController()));

  /*
   * Bundle potentials
   */
  auto pbePot = std::shared_ptr<PotentialBundle<SCFMode> >(
      new PBEPotentials<SCFMode>(
          activeSysPot,
          esiPot,
          naddXCPot,
          kin,
          ecpInt));

  // non-additive Dispersion correction
  auto nadDispCorrection = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
      supersystem->getGeometry(),settings.naddXCFunc);
  nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
      _environmentSystem->getGeometry(),settings.naddXCFunc);
  nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(settings.dispersion,
      _activeSystem->getGeometry(),settings.naddXCFunc);
  auto actSysDisp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(actsettings.dft.dispersion,
      _activeSystem->getGeometry(),_activeSystem->getSettings().dft.functional);

  /*
   * Calculate the frozen environment energy contributions and add them.
   */
  auto envEnergies = _environmentSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController();

  // add perturbative energy correction to the EnergyController.
  // TODO: If local MP2 is available, this may be different from zero.
  //       But will have to be calculate after the SCF for the active system
  //       in order to be able to project the occupied active orbitals.
  envEnergies->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION,0.0));

  std::shared_ptr<PotentialBundle<SCFMode> > envPot;
  double frozenEnvEnergy = 0.0;
  auto envsettings = _environmentSystem->getSettings();
  if (envsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
    envPot = _environmentSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
    envPot->getFockMatrix(densityMatrixEnvironment,envEnergies);
    frozenEnvEnergy = envEnergies->getEnergyComponent(ENERGY_CONTRIBUTIONS::HF_ENERGY);
  } else if (envsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    // Dispersion correction
    auto envDispEnergy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(_environmentSystem->getSettings().dft.dispersion,
        _environmentSystem->getGeometry(),_environmentSystem->getSettings().dft.functional);
    envEnergies->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION,envDispEnergy);
    envPot = _environmentSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
    envPot->getFockMatrix(densityMatrixEnvironment,envEnergies);
    frozenEnvEnergy = envEnergies->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
  } else {
    std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
    assert(false);
  }



//  const double frozenEnvEnergy = envEnergies->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);

  // Everything is set in environment ElectronicStructure -> save to file
  _environmentSystem->template getElectronicStructure<SCFMode>()->toHDF5(
      _environmentSystem->getHDF5BaseName(),_environmentSystem->getSettings().identifier);

  auto es = _activeSystem->getElectronicStructure<SCFMode>();
  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION,actSysDisp);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP,nadDispCorrection);
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_FROZEN_SUBSYSTEM_ENERGIES,frozenEnvEnergy);

  /*
   * Attach everything and run the SCF procedure
   */
  // manage output
  auto tmp1 = iOOptions.printFinalOrbitalEnergies;
  auto tmp2 = iOOptions.gridAccuracyCheck;
  iOOptions.printFinalOrbitalEnergies = false;
  iOOptions.gridAccuracyCheck = false;
  Scf<SCFMode>::perform(actsettings,es,pbePot);
  iOOptions.printFinalOrbitalEnergies = tmp1;
  iOOptions.gridAccuracyCheck = tmp2;

  //check for double hybrid functional
  auto functional = FunctionalClassResolver::resolveFunctional(actsettings.dft.functional);
  double MP2Correlation = 0.0;
  if (functional.isDoubleHybrid()){
    WarningTracker::printWarning("Warning: When using an double hybrid functional in embedding, the occupied environment orbitals are not explicitly projected out of the MP2 calculation! This may lead to errors.",iOOptions.printSCFCycleInfo);
    //perform MP2 for double hybrids
    if (actsettings.dft.densityFitting == Options::DENS_FITS::RI){
      _activeSystem->setBasisController(dynamic_pointer_cast<AtomCenteredBasisController>(
          supersystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL)),Options::BASIS_PURPOSES::AUX_CORREL);
      RIMP2<SCFMode> rimp2(_activeSystem, functional.getssScaling(), functional.getosScaling());
      MP2Correlation = rimp2.calculateCorrection();
    } else {
      assert(SCFMode==RESTRICTED && "MP2 is not available for unrestricted systems please use RI-MP2.");
      MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_activeSystem, functional.getssScaling(), functional.getosScaling());
      MP2Correlation = mp2EnergyCorrector.calculateElectronicEnergy();
    }
    MP2Correlation *= functional.getHfCorrelRatio();
  }

  // add energy to the EnergyController
  eCont->addOrReplaceComponent(
      std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION,MP2Correlation));

  if (iOOptions.printSCFResults){
    printSubSectionTitle("Final SCF Results");
    eCont->printAllComponents();
    if (iOOptions.printFinalOrbitalEnergies){
      print("");
      printSmallCaption("Orbital Energies");
      es->printMOEnergies();
      print("");
    }
  }
  return;
}

template<Options::SCF_MODES SCFMode>
void TDEmbeddingTask<SCFMode>::setActiveSystemBasis(){
  /*
   * For more information about the available truncation schemes see misc/BasisSetTruncationAlgorithms.h
   */
  _activeSystem->getGeometry()->printToFile(_activeSystem->getHDF5BaseName(),_activeSystem->getSettings().identifier);
  if (settings.truncAlgorithm == Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE) {
    /* =========================
     *   Use supersystem basis
     * =========================*/
    _activeSystem->getGeometry()->addAsDummy(*(_environmentSystem->getGeometry()));
    _activeSystem->getGeometry()->deleteIdenticalAtoms();
    /*
     * Reset bases, they will be rebuild when needed (using supersystem geometry)
     */
    _activeSystem->setBasisController(nullptr);
    _activeSystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_COULOMB);
    _activeSystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_CORREL);
  } else {
    BasisSetTruncationTask<SCFMode> basisSetTask(_activeSystem);
    basisSetTask.settings.truncationFactor = settings.truncationFactor;
    basisSetTask.settings.netThreshold = settings.netThreshold;
    basisSetTask.settings.truncAlgorithm = settings.truncAlgorithm;
    basisSetTask.run();
  }
}

template<Options::SCF_MODES SCFMode>
inline void TDEmbeddingTask<SCFMode>::checkInput() {
  if(settings.carterCycles!=0 and !settings.noSupRec){
    throw SerenityError(
        (string)"Zhang-Carter reconstruction does not support the double reconstruction feature!");
  }
  if(settings.useEnvSys and !settings.enforceCharges){
    throw SerenityError(
        (string)"Please enforce charges if you want to relax with respect to a precalculated environment!");
  }
  if(settings.useEnvSys and settings.truncAlgorithm != Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE){
    throw SerenityError(
        (string)"Relaxation with respect to precalc. environment is not compatible with basis truncation!");
  }
  //check for double hybrid functional
  auto envFunctional = FunctionalClassResolver::resolveFunctional(_environmentSystem->getSettings().dft.functional);
  if(envFunctional.isDoubleHybrid()) {
    throw SerenityError(
        (string)"Double hybrid functionals in the environment are not supported!");
  }
}
template<Options::SCF_MODES SCFMode>
inline std::shared_ptr<SystemController> TDEmbeddingTask<SCFMode>::setUpSupersystem() {
  /*
   * Load or create a supersystem from active and environment systems
   */
  std::shared_ptr<SystemController> supersystem(nullptr);
  std::shared_ptr<AtomCenteredGridController> supersystemGrid(nullptr);

  std::vector<std::shared_ptr<Atom> > superSystemAtoms;
  superSystemAtoms.insert(
      superSystemAtoms.end(), _activeSystem->getGeometry()->getAtoms().begin(), _activeSystem->getGeometry()->getAtoms().end());
  superSystemAtoms.insert(
      superSystemAtoms.end(), _environmentSystem->getGeometry()->getAtoms().begin(), _environmentSystem->getGeometry()->getAtoms().end());
  for (auto atom : _activeSystem->getGeometry()->getAtoms()){_activeAtoms.push_back(true);}
  for (auto atom : _environmentSystem->getGeometry()->getAtoms()){_activeAtoms.push_back(false);}
  int superSystemSpin = 0;
  int superSystemCharge = 0;

  if (settings.load == "" or settings.name == ""){
    superSystemSpin += _activeSystem->getSpin();
    superSystemCharge += _activeSystem->getCharge();
    superSystemSpin += _environmentSystem->getSpin();
    superSystemCharge += _environmentSystem->getCharge();


    // geometry of the entire system
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtoms);
    superSystemGeometry->deleteIdenticalAtoms();


    // supersystem grid TODO please check whether it is ok to use the active system settings and the default grid
    supersystemGrid = GridControllerFactory::produce(
        superSystemGeometry, _activeSystem->getSettings(), Options::GRID_PURPOSES::DEFAULT);

    /*
     * This takes the configuration of the environment system for the supersystem.
     * This approach is chosen for simplicity.
     */
    auto superSysSettings = _environmentSystem->getSettings();
    superSysSettings.spin = superSystemSpin;
    superSysSettings.charge = superSystemCharge;
    superSysSettings.name = _activeSystem->getSettings().name +"+"+superSysSettings.name;
    supersystem = std::make_shared<SystemController> (
        superSystemGeometry,
        superSysSettings);
    auto supSysBasisController=AtomCenteredBasisControllerFactory::produce(
        superSystemGeometry,
        supersystem->getSettings().basis.basisLibPath,
        supersystem->getSettings().basis.makeSphericalBasis,
        false,
        supersystem->getSettings().basis.firstECP,
        supersystem->getSettings().basis.label);

    supSysBasisController->getBasis();

    supersystem->setBasisController(supSysBasisController);

    /*
     * Supersystem SCF calculation
     */
    ScfTask<SCFMode> superSystemSCF(supersystem);
    superSystemSCF.run();
  }else{
    print((string)"Loading system: " + settings.name + " as supersystem");
    Settings superSysSettings;
    superSysSettings.load = settings.load;
    superSysSettings.name = settings.name;
    supersystem = std::make_shared<SystemController> (superSysSettings);
    superSystemSpin = supersystem->getSpin();
    superSystemCharge = supersystem->getCharge();
    supersystem->getAtomCenteredBasisController()->toHDF5(supersystem->getHDF5BaseName(),supersystem->getSettings().identifier);
    // Resort the active atoms vector.
    _activeAtoms = std::vector<bool>(supersystem->getGeometry()->getNAtoms(),false);
    for(const auto& actAtom : _activeSystem->getGeometry()->getAtoms()) {
      unsigned int indexInSupersystem = SystemSplittingTools<SCFMode>::matchAtom(supersystem->getGeometry(),actAtom);
      assert(indexInSupersystem < supersystem->getGeometry()->getNAtoms()&&"Error when loading supersystem. Active and supersystem atoms do not match!");
      _activeAtoms[indexInSupersystem] = true;
    }
    supersystemGrid =  GridControllerFactory::produce(
        supersystem->getGeometry(), _activeSystem->getSettings(), Options::GRID_PURPOSES::DEFAULT);
  }

  _activeSystem->setGridController(supersystemGrid);

  if(!settings.useEnvSys){

    _environmentSystem->getGeometry()->addAsDummy(*(_activeSystem->getGeometry()), true);
    _environmentSystem->getGeometry()->printToFile(_environmentSystem->getHDF5BaseName(),_environmentSystem->getSettings().identifier);
    _environmentSystem->getGeometry()->deleteIdenticalAtoms();

    auto envSysBasisController=AtomCenteredBasisControllerFactory::produce(
        _environmentSystem->getGeometry(),
        _environmentSystem->getSettings().basis.basisLibPath,
        _environmentSystem->getSettings().basis.makeSphericalBasis,
        false,
        _environmentSystem->getSettings().basis.firstECP,
        "");
    envSysBasisController->toHDF5(_environmentSystem->getHDF5BaseName(),_environmentSystem->getSettings().identifier);

    _environmentSystem->setGridController(supersystemGrid);
    _environmentSystem->setBasisController(envSysBasisController);
    _environmentSystem->setBasisController(nullptr,
        Options::BASIS_PURPOSES::AUX_COULOMB);
  }else{
    supersystem->setBasisController(dynamic_pointer_cast<AtomCenteredBasisController>(
        _environmentSystem->getBasisController()));
    supersystem->setBasisController(dynamic_pointer_cast<AtomCenteredBasisController>(
        _environmentSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
        Options::BASIS_PURPOSES::AUX_COULOMB);
    _activeAtoms.clear();
    for (auto atom : _environmentSystem->getGeometry()->getAtoms()){if(!atom->isDummy())_activeAtoms.push_back(false);}
    for (auto atom : _activeSystem->getGeometry()->getAtoms()){if(!atom->isDummy())_activeAtoms.push_back(true);}
  }
  return supersystem;
}

template<>
inline SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::string>
TDEmbeddingTask<Options::SCF_MODES::RESTRICTED>::getDensityNameSuffix() {
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, std::string> suffix("");
  return suffix;
}

template<>
inline SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::string>
TDEmbeddingTask<Options::SCF_MODES::UNRESTRICTED>::getDensityNameSuffix() {
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, std::string> suffix("");
  suffix.alpha = "alpha";
  suffix.beta = "beta";
  return suffix;
}


template class TDEmbeddingTask<Options::SCF_MODES::RESTRICTED>;
template class TDEmbeddingTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
