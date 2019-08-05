/**
 * @file FaTConvergenceAccelerator.cpp
 *
 * @date Mar 29, 2018
 * @author Moritz Bensberg
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
#include "misc/FaTConvergenceAccelerator.h"

/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "grid/GridControllerFactory.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/diis/DIIS.h"
#include "misc/Timing.h"
#include "memory/MemoryManager.h"
#include "potentials/BUReconstructionPotential.h"

/* Construction of embedded fock matrix. */
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/HuzinagaFDEProjectionPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/ECPInteractionPotential.h"
#include "potentials/LevelshiftHybridPotential.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FaTConvergenceAccelerator<SCFMode>::FaTConvergenceAccelerator(
    unsigned int maxStore,
    double conditionNumberThreshold,
    const FreezeAndThawTaskSettings& settings,
    std::vector<std::shared_ptr<SystemController> > activeSystems,
    std::vector<std::shared_ptr<SystemController> > environmentSystems):
    _activeSystems(activeSystems),
    _environmentSystems(environmentSystems),
    _settings(settings){
  // initialize the error vector, diis and density vector
  // calculate cache limit for the vector controller with a buffer of 2 GB
  auto memoManger = MemoryManager::getInstance();
  long long availableMemory =
      (memoManger->getAvailableSystemMemory()>2e+9)?memoManger->getAvailableSystemMemory()-2e+9 : 0;
  double cacheLimit = (double) availableMemory/((double)maxStore+3.0);
  unsigned int vectorLength = 0;
  for (unsigned int i = 0; i < _activeSystems.size(); ++i) {
    unsigned int nBasisFunctionsI = _activeSystems[i]->getBasisController()->getNBasisFunctions();
    vectorLength += nBasisFunctionsI*nBasisFunctionsI;
  }
  _densityVector = std::make_shared<VectorOnDiskStorageController>(cacheLimit,"Density.h5");
  _fpsMinusSPF = std::make_shared<VectorOnDiskStorageController>(cacheLimit,"Error.h5");
  _diis = std::make_shared<DIIS>(maxStore, conditionNumberThreshold, false, true);
}
template<Options::SCF_MODES SCFMode>
void FaTConvergenceAccelerator<SCFMode>::accelerateConvergence(double energy) {
  takeTime("Freeze and Thaw DIIS");
  _cycle++;
  // map the density matrices to a vector
  for (auto& sys : _activeSystems) {
    unsigned int nBasisFunctions = sys->getBasisController()->getNBasisFunctions();
    DensityMatrix<SCFMode> p = sys->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix();
    unsigned int vectorSizeFactor = (SCFMode == Options::SCF_MODES::RESTRICTED)? 1 : 2;
    auto newVectorSegment = std::make_shared<Eigen::VectorXd> (vectorSizeFactor*nBasisFunctions*nBasisFunctions);
    unsigned int lastSegmentIndex = 0;
    for_spin(p) {
      newVectorSegment->block(lastSegmentIndex,0,nBasisFunctions*nBasisFunctions,1)
                      = Eigen::Map<Eigen::VectorXd>(p_spin.data(),nBasisFunctions*nBasisFunctions);
      lastSegmentIndex += nBasisFunctions*nBasisFunctions;
    };
    _densityVector->storeVectorSegment(newVectorSegment,sys->getSettings().name);
  }
  // Calculate RMSD of the density
  double rmsd = calcRMSDofDensity();

  // Start DIIS if the RMSD is below a given threshold.
  // Else damp.
  std::cout << "-------------------------------------"<<std::endl;
  std::cout << "Density RMSD (all systems): " << rmsd << std::endl;
  if (rmsd < _settings.diisStart && rmsd > _settings.diisEnd ) {
    // Calculate error measure.
    calcFPSminusSPF();
    std::cout << "+++ performing DIIS step +++" << std::endl;
    std::cout << "-------------------------------------"<<std::endl;
    // perform DIIS step
    _diis->optimize(energy, *_densityVector, *_fpsMinusSPF);

    // Set the optimized density in the system controller
    for (auto& sys : _activeSystems) {
      DensityMatrix<SCFMode> optDensity(sys->getBasisController());
      unsigned int nBasisFunctions = sys->getBasisController()->getNBasisFunctions();
      auto optVectorSegmentD = _densityVector->getVectorSegment(sys->getSettings().name);
      unsigned int lastSegmentIndex = 0;
      for_spin(optDensity) {
        optDensity_spin = Eigen::Map<Eigen::MatrixXd>(
            optVectorSegmentD->block(lastSegmentIndex,0,nBasisFunctions*nBasisFunctions,1).data(),
            nBasisFunctions,nBasisFunctions);
        lastSegmentIndex+=nBasisFunctions*nBasisFunctions;
      };
      sys->getElectronicStructure<SCFMode>()->getDensityMatrixController()->setDensityMatrix(optDensity);

    } /* for sys */
  } /* rmsd < _diisStartError */
  else {
    std::cout << "-------------------------------------"<<std::endl;
  }
  timeTaken(3,"Freeze and Thaw DIIS");
}

template<Options::SCF_MODES SCFMode>
void FaTConvergenceAccelerator<SCFMode>::calcFPSminusSPF() {
  /*
   * 1. Loop over active systems.
   * 2. Calculate embedded fock matrix.
   * 3. Calculate FPS-SPF
   * 4. Map the matrix to an vector segment and append the block
   */
  for (unsigned int i = 0; i < _activeSystems.size(); ++i) {
    auto f = calcEmbeddedFockMatrix(i);
    auto p = _activeSystems[i]->getElectronicStructure<SCFMode>()->getDensityMatrix();
    auto s = _activeSystems[i]->getOneElectronIntegralController()->getOverlapIntegrals();
    unsigned int nBasisFunctions = _activeSystems[i]->getBasisController()->getNBasisFunctions();

    unsigned int vectorSizeFactor = (SCFMode == Options::SCF_MODES::RESTRICTED)? 1 : 2;
    auto newVectorSegmentFPS_SPF = std::make_shared<Eigen::VectorXd>(vectorSizeFactor*nBasisFunctions*nBasisFunctions);
    auto newVectorSegmentF = std::make_shared<Eigen::VectorXd>(vectorSizeFactor*nBasisFunctions*nBasisFunctions);
    unsigned int lastSegmentIndex = 0;
    for_spin(f,p) {
      FockMatrix<Options::SCF_MODES::RESTRICTED> fps_spf(_activeSystems[i]->getBasisController());
      fps_spf = f_spin * p_spin * s - s * p_spin * f_spin;
      newVectorSegmentFPS_SPF->block(lastSegmentIndex,0,s.rows()*s.cols(),1)
                             = Eigen::Map<Eigen::VectorXd>(fps_spf.data(),s.rows()*s.cols());
      // Save the fock matrix
      lastSegmentIndex += s.rows()*s.cols();
    };
    _fpsMinusSPF->storeVectorSegment(newVectorSegmentFPS_SPF,_activeSystems[i]->getSettings().name);
  }
}
template<Options::SCF_MODES SCFMode>
double FaTConvergenceAccelerator<SCFMode>::calcRMSDofDensity() {
  double rmsd = 0.0;
  // If an old density vector is available calculate the RMSD.
  // Otherwise return inf.
  if (_oldDensityVector) {
    double sumOfSquares = 0.0;
    // Loop over systems
    for (const auto& sys : _activeSystems) {
      std::string label = sys->getSettings().name;
      auto differenceSegment =
          *_densityVector->getVectorSegment(label)
          -*_oldDensityVector->getVectorSegment(label);
      double tmp = differenceSegment.cwiseProduct(differenceSegment).sum();
      sumOfSquares += tmp;
    }
    rmsd = std::sqrt(sumOfSquares/_densityVector->size());
    // Set new "old" density vector.
    _oldDensityVector = nullptr;
    _oldDensityVector.reset(new VectorOnDiskStorageController (*_densityVector,"OldDensity.h5"));
  } else {
    rmsd = std::numeric_limits<double>::infinity();
    _oldDensityVector = std::make_shared<VectorOnDiskStorageController> (*_densityVector,"OldDensity.h5");
  }
  return rmsd;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> FaTConvergenceAccelerator<SCFMode>::calcEmbeddedFockMatrix(
    unsigned int activeSystemIndex) {
  // list of environment systems for this active system.
  std::vector<std::shared_ptr<SystemController> > envSystems;
  for (unsigned int i = 0; i < _activeSystems.size(); ++i) {
    if (i != activeSystemIndex) envSystems.push_back(_activeSystems[i]);
  }
  for (unsigned int i = 0 ; i < _environmentSystems.size(); ++i) {
    envSystems.push_back(_environmentSystems[i]);
  }
  // list of environment density matrices (their controllers)
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDensities;
  for (auto& sys : envSystems){
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
  } /* for i */
  /* Core part of the fock matrix */
  auto activeSystem = _activeSystems[activeSystemIndex];
  auto actsettings = activeSystem->getSettings();
  std::shared_ptr<PotentialBundle<SCFMode> > activeSysPot;
  if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
    activeSysPot = activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  } else if (actsettings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    activeSysPot = activeSystem->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>();
  } else {
    std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
    assert(false);
  }
  // ECP Contribution
  std::vector<std::shared_ptr<Atom> > envAtoms;
  for (const auto& sys : envSystems) {
    for (const auto& atom : sys->getGeometry()->getAtoms()) {
      envAtoms.push_back(atom);
    }
  }
  std::shared_ptr<Potential<SCFMode> > ecpInt
  (new ECPInteractionPotential<SCFMode>(
      activeSystem,
      activeSystem->getGeometry()->getAtoms(),
      envAtoms,
      envDensities,
      activeSystem->getBasisController()));
  /* Electrostatic embedding */
  // geometries of the environment subsystems
  std::vector<std::shared_ptr<const Geometry> > envGeometries;
  for (auto sys : envSystems){
    envGeometries.push_back(sys->getGeometry());
  }
  std::shared_ptr<PotentialBundle<SCFMode> > esiPot;
  if (actsettings.dft.densityFitting==Options::DENS_FITS::RI){
    std::vector<std::shared_ptr<BasisController> > envAuxBasis;
    for(auto& env : envSystems){
      envAuxBasis.push_back(env->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
    }
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            activeSystem,
            envSystems,
            activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            activeSystem->getGeometry(),
            envDensities,
            envGeometries,
            activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
            envAuxBasis));
  } else {
    esiPot = std::shared_ptr<PotentialBundle<SCFMode> >(
        new ESIPotentials<SCFMode>(
            activeSystem,
            envSystems,
            activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
            activeSystem->getGeometry(),
            envDensities,
            envGeometries));
  }
  auto es = activeSystem->getElectronicStructure<SCFMode>();
  auto grid = activeSystem->getGridController();
  auto eCont = es->getEnergyComponentController();
  eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_EXACT_EXCHANGE,0.0);

  std::shared_ptr<NAddFuncPotential<SCFMode> > naddKinPot;
  auto naddXCfunc = FunctionalClassResolver::resolveFunctional(_settings.naddXCFunc);
  auto naddXCPot = std::make_shared<NAddFuncPotential<SCFMode> >(
      activeSystem,
      activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      envDensities,
      grid,
      naddXCfunc,
      std::make_pair(true,std::vector<std::shared_ptr<EnergyComponentController> >(0)),
      false,
      false);

  std::shared_ptr<Potential<SCFMode> > kin;
  if(_settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NADD_FUNC){
    kin = std::make_shared<NAddFuncPotential<SCFMode> >(
        activeSystem,
        activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        grid,
        FunctionalClassResolver::resolveFunctional(_settings.naddKinFunc),
        std::make_pair(true,std::vector<std::shared_ptr<EnergyComponentController> >(0)),
        false,
        false);
  }else if(_settings.embeddingMode==Options::KIN_EMBEDDING_MODES::LEVELSHIFT){
    kin = std::make_shared<LevelshiftHybridPotential<SCFMode> > (
        activeSystem,
        envSystems,
        _settings.levelShiftParameter,
        (_settings.basisFunctionRatio == 0.0)? false : true,
         _settings.basisFunctionRatio,
         _settings.borderAtomThreshold,
         grid,
         FunctionalClassResolver::resolveFunctional(_settings.naddKinFunc),
         false);
  }else if(_settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HUZINAGA){
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode> >(
        activeSystem,
        envSystems,
        naddXCfunc,
        _settings.truncateProjector,
        _settings.projecTruncThresh,
        false,
        false,
        _settings.distantKinFunc,
        _settings.basisFunctionRatio,
        _settings.borderAtomThreshold,
        grid,
        _settings.naddKinFunc,
        0,//_settings.gridCutOff,
        std::vector<std::shared_ptr<EnergyComponentController> >(0));
  }else if(_settings.embeddingMode==Options::KIN_EMBEDDING_MODES::HOFFMANN){
    kin = std::make_shared<HuzinagaFDEProjectionPotential<SCFMode> >(
        activeSystem,
        envSystems,
        naddXCfunc,
        _settings.truncateProjector,
        _settings.projecTruncThresh,
        true,
        false,
        _settings.distantKinFunc,
        _settings.basisFunctionRatio,
        _settings.borderAtomThreshold,
        grid,
        _settings.naddKinFunc,
        _settings.gridCutOff,
        std::vector<std::shared_ptr<EnergyComponentController> >(0));
  }else if(_settings.embeddingMode==Options::KIN_EMBEDDING_MODES::RECONSTRUCTION){
    kin=std::make_shared<BUReconstructionPotential<SCFMode>>(
        activeSystem,
        activeSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        envDensities,
        grid,
        naddXCfunc,
        _environmentSystems,
        _settings.smoothFactor,
        _settings.potentialBasis.empty()? activeSystem->getBasisController()->getBasisString() : _settings.potentialBasis,
        _settings.singValThreshold,
        _settings.lbDamping,
        _settings.lbCycles,
        _settings.carterCycles);
  }else{
    kin = std::shared_ptr<Potential<SCFMode> >(
        new ZeroPotential<SCFMode>(activeSystem->getBasisController()));
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

  return fdePot->getFockMatrix(activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix(),
      std::make_shared<EnergyComponentController>());
}

template class FaTConvergenceAccelerator<Options::SCF_MODES::RESTRICTED>;
template class FaTConvergenceAccelerator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
