/**
 * @file FreezeAndThawTask.cpp
 *
 * @date Oct 15, 2015
 * @author Jan Unsleber
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
#include "tasks/FreezeAndThawTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/BasisExtension.h"
#include "math/optimizer/BFGS.h"
#include "data/grid/CustomDensityOnGridController.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/DensityOnGridFactory.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "tasks/FDETask.h"
#include "io/FormattedOutput.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "tasks/GradientTask.h"
#include "grid/GridControllerFactory.h"
#include "tasks/HessianTask.h"
#include "io/IOOptions.h"
#include "integrals/wrappers/Libint.h"
#include "math/optimizer/Optimizer.h"
#include "tasks/ScfTask.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "misc/Timing.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "misc/FaTConvergenceAccelerator.h"
/* Include Std and External Headers */
#include <stdio.h>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FreezeAndThawTask<SCFMode>::FreezeAndThawTask(
    const std::vector<std::shared_ptr<SystemController> >& activeSystems,
    const std::vector<std::shared_ptr<SystemController> >& passiveSystems) :
       _activeSystems(activeSystems),
       _passiveSystems(passiveSystems){
}

template<Options::SCF_MODES SCFMode>
void FreezeAndThawTask<SCFMode>::run() {
  if (settings.extendBasis) {
    // increase the basis for all active systems
	BasisExtension basisExtender;
	basisExtender.extendAllSystems(_activeSystems, settings.basisExtThresh);
	std::cout << "Extended all basis sets!" << std::endl<<std::endl;
  }
  // shorten the output
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool gridAcc(iOOptions.gridAccuracyCheck);
  int timings(iOOptions.timingsPrintLevel);
  bool printTitle = true;

  switch (settings.printLevel) {
    case 0:
      iOOptions.printSCFCycleInfo = false;
      iOOptions.printSCFResults = false;
      iOOptions.printGridInfo = false;
      printTitle = false;
      break;
    case 1:
      iOOptions.printSCFCycleInfo = false;
      iOOptions.printSCFResults = true;
      break;
    case 2:
      iOOptions.printSCFCycleInfo = true;
      iOOptions.printSCFResults = true;
      break;
  }
  iOOptions.gridAccuracyCheck = false;
  iOOptions.timingsPrintLevel = 1;

  // keep Libint throughout freeze and thaw cycles
  auto& libint = Libint::getInstance();
  libint.keepEngines(libint2::Operator::kinetic,0,2);
  libint.keepEngines(libint2::Operator::nuclear,0,2);
  libint.keepEngines(libint2::Operator::overlap,0,2);
  libint.keepEngines(libint2::Operator::coulomb,0,2);
  libint.keepEngines(libint2::Operator::coulomb,0,3);
  libint.keepEngines(libint2::Operator::coulomb,0,4);

  if(settings.makeSuperSystemBasis){
    /*
     * Reset bases to whatever is requested in userinput
     * to store basis information in atoms
     */
    for(auto system : _activeSystems){
      system->setBasisController(nullptr);
      system->getBasisController()->getBasis();
    }
    for(auto system : _passiveSystems){
      system->setBasisController(nullptr);
      system->getBasisController()->getBasis();
    }
    for(auto system1 : _activeSystems){
      for(auto system2 : _activeSystems){
        if(system1==system2) continue;
        system1->getGeometry()->addAsDummy(*system2->getGeometry());
      }
      for(auto system2 : _passiveSystems){
        if(system1==system2) continue;
        system1->getGeometry()->addAsDummy(*system2->getGeometry());
      }
    }
    for(auto system1 : _passiveSystems){
      for(auto system2 : _activeSystems){
        if(system1==system2) continue;
        system1->getGeometry()->addAsDummy(*system2->getGeometry());
      }
      for(auto system2 : _passiveSystems){
        if(system1==system2) continue;
        system1->getGeometry()->addAsDummy(*system2->getGeometry());
      }
    }
    for(auto system : _activeSystems){
      auto sysBasisController=AtomCenteredBasisControllerFactory::produce(
          system->getGeometry(),
          system->getSettings().basis.basisLibPath,
          system->getSettings().basis.makeSphericalBasis,
          false,
          system->getSettings().basis.firstECP,
          "");
      system->setBasisController(sysBasisController);
    }
    for(auto system : _passiveSystems){
      auto sysBasisController=AtomCenteredBasisControllerFactory::produce(
          system->getGeometry(),
          system->getSettings().basis.basisLibPath,
          system->getSettings().basis.makeSphericalBasis,
          false,
          system->getSettings().basis.firstECP,
          "");
      system->setBasisController(sysBasisController);
    }
  }

  // Isolated runs if needed
  for (auto sys : _activeSystems){
    if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
      if (printTitle && !sys->template hasElectronicStructure<Options::SCF_MODES::RESTRICTED>()){
        printSubSectionTitle((std::string)"Isolated Run");
      }
      sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
    } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED){
      if (printTitle && !sys->template hasElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()){
        printSubSectionTitle((std::string)"Isolated Run");
      }
      sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
    } else {
      assert(false);
    }
  }
  for (auto sys : _passiveSystems){
    if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
      sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
    } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED){
      sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
    } else {
      assert(false);
    }
  }

  double superSysE = 0.0;

  // Supersystem grid
  std::shared_ptr<GridController> supersystemgrid(nullptr);
  std::vector<std::shared_ptr<DensityOnGridController<SCFMode> > >densOnGridControllers(
      _activeSystems.size()+_passiveSystems.size(),nullptr);

  std::shared_ptr<FaTConvergenceAccelerator<SCFMode> > fatConvAcc = nullptr;
  if(settings.useConvAcceleration) fatConvAcc = std::make_shared<FaTConvergenceAccelerator<SCFMode> >(
      _activeSystems[0]->getSettings().scf.diisMaxStore,
      _activeSystems[0]->getSettings().scf.diisConditionNumberThreshold,
      settings,
      _activeSystems,
      _passiveSystems);

  // Freeze and thaw cycles
  for (unsigned int cycle=0; cycle<settings.maxCycles;cycle++ ){

    // Ouput and init
    if (printTitle) printSubSectionTitle((std::string)"Freeze-and-Thaw Cycle: "+(cycle+1));
    std::vector<bool> convergedSubsystems(_activeSystems.size(),false);

    // Active system cycles
    for (unsigned int nSystem=0; nSystem<_activeSystems.size();nSystem++){

      // Output and init
      auto activeSystem =  _activeSystems[nSystem];

      // Keep densities on grid alive
      if (nSystem>0 and settings.gridCutOff<0.0){
        for (unsigned int nSystem=0; nSystem<_activeSystems.size();nSystem++){
          densOnGridControllers[nSystem] =
              DensityOnGridFactory<SCFMode>::produce(
                  _activeSystems[nSystem]->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
                  supersystemgrid,1,activeSystem->getSettings());
        }
        for (unsigned int nSystem=0; nSystem<_passiveSystems.size();nSystem++){
          densOnGridControllers[nSystem+_activeSystems.size()] =
              DensityOnGridFactory<SCFMode>::produce(
                  _passiveSystems[nSystem]->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
                  supersystemgrid,1,activeSystem->getSettings());
        }
      }

      if (printTitle) printBigCaption((std::string)"Active System: "+activeSystem->getSettings().name);
      // Store old energy of the active System
      DensityMatrix<RESTRICTED> rhoOld(activeSystem->getBasisController());
      if (activeSystem->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        rhoOld = activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix().total();
      } else if (activeSystem->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        rhoOld = activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
      } else {
        assert(false);
      }
      // Set up systems
      auto passiveSystems = _activeSystems;
      passiveSystems.erase(std::remove(passiveSystems.begin(), passiveSystems.end(), activeSystem), passiveSystems.end());
      passiveSystems.insert(passiveSystems.end(), _passiveSystems.begin(), _passiveSystems.end());
      // Update the active System
      if (activeSystem->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        FDETask<Options::SCF_MODES::RESTRICTED> task(activeSystem,passiveSystems);
        task.settings.naddKinFunc = settings.naddKinFunc;
        task.settings.naddXCFunc = settings.naddXCFunc;
        task.settings.embeddingMode = settings.embeddingMode;
        task.settings.gridCutOff = settings.gridCutOff;
        task.settings.dispersion = settings.dispersion;
        task.settings.smoothFactor = settings.smoothFactor;
        task.settings.potentialBasis = settings.potentialBasis;
        task.settings.singValThreshold = settings.singValThreshold;
        task.settings.lbDamping = settings.lbDamping;
        task.settings.lbCycles = settings.lbCycles;
        task.settings.carterCycles = settings.carterCycles;
        task.settings.smallSupersystemGrid = settings.smallSupersystemGrid;
        task.settings.finalGrid=false;
        task.settings.truncateProjector = settings.truncateProjector;
        task.settings.projecTruncThresh = settings.projecTruncThresh;
        task.settings.distantKinFunc = settings.distantKinFunc;
        if (settings.gridCutOff<0) task.setSuperSystemGrid(supersystemgrid);
        task.run();
        supersystemgrid = task.getSuperSystemGrid();
      } else if (activeSystem->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        FDETask<Options::SCF_MODES::UNRESTRICTED> task(activeSystem,passiveSystems);
        task.settings.naddKinFunc = settings.naddKinFunc;
        task.settings.naddXCFunc = settings.naddXCFunc;
        task.settings.embeddingMode = settings.embeddingMode;
        task.settings.gridCutOff = settings.gridCutOff;
        task.settings.dispersion = settings.dispersion;
        task.settings.smoothFactor = settings.smoothFactor;
        task.settings.potentialBasis = settings.potentialBasis;
        task.settings.singValThreshold = settings.singValThreshold;
        task.settings.lbDamping = settings.lbDamping;
        task.settings.lbCycles = settings.lbCycles;
        task.settings.carterCycles = settings.carterCycles;
        task.settings.smallSupersystemGrid = settings.smallSupersystemGrid;
        task.settings.finalGrid=false;
        task.settings.truncateProjector = settings.truncateProjector;
        task.settings.projecTruncThresh = settings.projecTruncThresh;
        task.settings.distantKinFunc = settings.distantKinFunc;
        if (settings.gridCutOff<0) task.setSuperSystemGrid(supersystemgrid);
        task.run();
        supersystemgrid = task.getSuperSystemGrid();
      } else {
        assert(false);
      }
      // Converged subsystem?
      double deltaRho = std::numeric_limits<double>::infinity();
      if (activeSystem->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        deltaRho = (rhoOld - activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix().total()).array().square().sum();
        deltaRho = sqrt(deltaRho);
      } else if (activeSystem->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        deltaRho = (rhoOld - activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total()).array().square().sum();
        deltaRho = sqrt(deltaRho);
      } else {
        assert(false);
      }
      if (fabs(deltaRho) < settings.convThresh){
        convergedSubsystems[nSystem] = true;
      }

      // Calculate total energies
      assert((activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
    		  activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF )&&
			  "Unknown electronic theory!");
      auto energyContribution = (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)?
    		  ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT : ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT;
      if (activeSystem->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        superSysE =activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energyContribution);
      } else if (activeSystem->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        superSysE =activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energyContribution);
      } else {
        assert(false);
      }
      // Output
      if (printTitle) printf("\n");
      if (printTitle) printSmallCaption("Freeze-and-Thaw Info");
      if (printTitle) printf("Freeze-and-Thaw Cycle               %16d\n",cycle+1);
      if (printTitle) printf("Active subsystem:                   %16d\n",nSystem+1);
      if (printTitle) printf("Density change (abs. deviation):        %e\n",deltaRho);
      if (printTitle) printf("Total supersystem energy:         %18.10f\n",superSysE);
      if (printTitle) printf("\n");
    }

    // Convergence check.
    if (std::find(convergedSubsystems.begin(), convergedSubsystems.end(), false) == convergedSubsystems.end()){
      printSubSectionTitle((std::string)"Converged after "+(cycle+1)+" freeze-and-thaw cycles!");
      break;
    }
    // Accelerate convergence
    if(settings.useConvAcceleration) fatConvAcc->accelerateConvergence(superSysE);
  }


  /* =======================================
   *   Final Energy Evaluation (if needed)
   * ======================================= */
  if (settings.gridCutOff>=0.0 and settings.embeddingMode==Options::KIN_EMBEDDING_MODES::NADD_FUNC){
    printSmallCaption("Creating final grid for energy evaluation");
    // atoms of all subsystems
    auto finalGridGeometry=std::make_shared<Geometry>();

    for (auto sys : _activeSystems){
      *finalGridGeometry+=*sys->getGeometry();
    }

    for (auto sys : _passiveSystems){
      *finalGridGeometry+=*sys->getGeometry();
    }
    finalGridGeometry->deleteIdenticalAtoms();

    Timings::takeTime("FaT - Final Energy Eval.");
    Options::GRID_PURPOSES finalGridacc = (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL:Options::GRID_PURPOSES::DEFAULT;
    auto finalGrid = GridControllerFactory::produce(
        finalGridGeometry, _activeSystems[0]->getSettings(), finalGridacc);

    auto xcfunc = FunctionalClassResolver::resolveFunctional(settings.naddXCFunc);
    auto kinefunc = FunctionalClassResolver::resolveFunctional(settings.naddKinFunc);
    bool gga(kinefunc.getFunctionalClass()!=FUNCTIONAL_CLASSES::LDA or
              xcfunc.getFunctionalClass()!=FUNCTIONAL_CLASSES::LDA);

    std::unique_ptr<DensityOnGrid<SCFMode> > superSysDens(new DensityOnGrid<SCFMode>(finalGrid));
    std::unique_ptr<Gradient<DensityOnGrid<SCFMode> > > superSysdensGrad(nullptr);
    if (gga) superSysdensGrad = std::move(makeGradientPtr<DensityOnGrid<SCFMode> >(finalGrid));

    double naddXcEnergy = 0.0;
    double naddKinEnergy = 0.0;

    for (unsigned int nSystem=0; nSystem<_activeSystems.size();nSystem++){
      auto calc = DensityOnGridFactory<SCFMode>::produce(
          _activeSystems[nSystem]->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
          finalGrid,1,_activeSystems[nSystem]->getSettings());
      *superSysDens       += calc->getDensityOnGrid();
      if (gga){
        superSysdensGrad->x += calc->getDensityGradientOnGrid().x;
        superSysdensGrad->y += calc->getDensityGradientOnGrid().y;
        superSysdensGrad->z += calc->getDensityGradientOnGrid().z;
      }

      XCFun<SCFMode> xcfun(128);
      naddXcEnergy -= xcfun.calcData(
          FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
          xcfunc,
          calc,0).energy;
      naddKinEnergy -= xcfun.calcData(
          FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
          kinefunc,
          calc,0).energy;
     }

    for (unsigned int nSystem=0; nSystem<_passiveSystems.size();nSystem++){
      auto calc = DensityOnGridFactory<SCFMode>::produce(
          _passiveSystems[nSystem]->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
          finalGrid,1,_passiveSystems[nSystem]->getSettings());
      *superSysDens       += calc->getDensityOnGrid();
      if (gga){
        superSysdensGrad->x += calc->getDensityGradientOnGrid().x;
        superSysdensGrad->y += calc->getDensityGradientOnGrid().y;
        superSysdensGrad->z += calc->getDensityGradientOnGrid().z;
      }
      XCFun<SCFMode> xcfun(128);
      naddXcEnergy -= xcfun.calcData(
          FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
          xcfunc,
          calc,0).energy;
      naddKinEnergy -= xcfun.calcData(
          FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
          kinefunc,
          calc,0).energy;
    }

    auto superSysDensOnGrid = std::make_shared<CustomDensityOnGridController<SCFMode> >(
        std::move(superSysDens),std::move(superSysdensGrad));
    XCFun<SCFMode> xcfun(128);
    naddXcEnergy += xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        xcfunc,
        superSysDensOnGrid,0).energy;
    naddKinEnergy += xcfun.calcData(
        FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS,
        kinefunc,
        superSysDensOnGrid,0).energy;

    for(auto& sys : _activeSystems){
      auto eCont=sys->getElectronicStructure<SCFMode>()->getEnergyComponentController();
      eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC,naddXcEnergy);
      eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC,naddKinEnergy);
    }
    Timings::timeTaken("FaT - Final Energy Eval.");
  }
  std::cout << std::endl;
  std::cout << std::endl;
  printSubSectionTitle("Final Freeze-and-Thaw Energies");
  for (unsigned int nSystem=0; nSystem<_activeSystems.size();nSystem++){
    printBigCaption((std::string)"Active System: "+_activeSystems[nSystem]->getSettings().name);
    auto eCont=_activeSystems[nSystem]->getElectronicStructure<SCFMode>()->getEnergyComponentController();
    eCont->printAllComponents();
    std::cout << std::endl;
    std::cout << std::endl;
  }

  // allow for freeing of Libint engines
  libint.freeEngines(libint2::Operator::kinetic,0,2);
  libint.freeEngines(libint2::Operator::nuclear,0,2);
  libint.freeEngines(libint2::Operator::overlap,0,2);
  libint.freeEngines(libint2::Operator::coulomb,0,2);
  libint.freeEngines(libint2::Operator::coulomb,0,3);
  libint.freeEngines(libint2::Operator::coulomb,0,4);

  // reset print options
  iOOptions.printSCFCycleInfo =info;
  iOOptions.printSCFResults =results;
  iOOptions.gridAccuracyCheck = gridAcc;
  iOOptions.timingsPrintLevel = timings;
}

template class FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>;
template class FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
