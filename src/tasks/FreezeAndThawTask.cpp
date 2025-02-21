/**
 * @file FreezeAndThawTask.cpp
 *
 * @date Oct 15, 2015
 * @author Jan Unsleber
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
#include "tasks/FreezeAndThawTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/BasisExtension.h"
#include "data/ElectronicStructure.h"
#include "data/grid/CustomDensityOnGridController.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/DensityOnGridFactory.h"
#include "data/grid/ElectrostaticPotentialOnGridController.h" //Potential on grid construction.
#include "data/grid/SupersystemDensityOnGridController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "geometry/MolecularSurfaceController.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "io/IOOptions.h"
#include "math/optimizer/BFGS.h"
#include "math/optimizer/Optimizer.h"
#include "misc/FaTConvergenceAccelerator.h"
#include "misc/Timing.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <stdio.h>
#include <string>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FreezeAndThawTask<SCFMode>::FreezeAndThawTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                                              const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void FreezeAndThawTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _activeSystems);
  if (settings.extendBasis) {
    // increase the basis for all active systems
    BasisExtension basisExtender;
    basisExtender.extendAllSystems(_activeSystems, settings.basisExtThresh);
    std::cout << "Extended all basis sets!" << std::endl << std::endl;
  }
  // shorten the output
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool gridAcc(iOOptions.gridAccuracyCheck);
  bool printTitle = true;

  switch (GLOBAL_PRINT_LEVEL) {
    case Options::GLOBAL_PRINT_LEVELS::MINIMUM:
      iOOptions.printSCFCycleInfo = false;
      iOOptions.printSCFResults = false;
      iOOptions.printGridInfo = false;
      printTitle = false;
      break;
    case Options::GLOBAL_PRINT_LEVELS::NORMAL:
    case Options::GLOBAL_PRINT_LEVELS::VERBOSE:
    case Options::GLOBAL_PRINT_LEVELS::DEBUGGING:
      iOOptions.printSCFCycleInfo = true;
      iOOptions.printSCFResults = true;
      break;
  }
  iOOptions.gridAccuracyCheck = false;

  // keep Libint throughout freeze and thaw cycles
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::kinetic, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::nuclear, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::overlap, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  // Isolated runs if needed
  for (auto sys : _activeSystems) {
    if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
      if (printTitle && !sys->template hasElectronicStructure<Options::SCF_MODES::RESTRICTED>()) {
        printSubSectionTitle((std::string) "Isolated Run");
      }
      sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
    }
    else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
      if (printTitle && !sys->template hasElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()) {
        printSubSectionTitle((std::string) "Isolated Run");
      }
      sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
    }
    else {
      assert(false);
    }
  }
  for (auto sys : _passiveSystems) {
    if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
      sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
    }
    else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
      sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
    }
    else {
      assert(false);
    }
  }

  double superSysE = 0.0;

  // Supersystem grid
  std::shared_ptr<GridController> supersystemgrid(nullptr);
  std::vector<std::shared_ptr<DensityOnGridController<SCFMode>>> densOnGridControllers(
      _activeSystems.size() + _passiveSystems.size(), nullptr);

  std::shared_ptr<FaTConvergenceAccelerator<SCFMode>> fatConvAcc = nullptr;
  if (settings.useConvAcceleration) {
    if (settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::ALMO)
      WarningTracker::printWarning("Warning: Use of convergence acceleration is not supported for ALMOs.",
                                   iOOptions.printSCFCycleInfo);
    else
      fatConvAcc = std::make_shared<FaTConvergenceAccelerator<SCFMode>>(_activeSystems[0]->getSettings().scf.diisMaxStore,
                                                                        settings, _activeSystems, _passiveSystems);
  }

  bool molecularSurfaceIsInitialized = false;
  // Freeze and thaw cycles
  for (unsigned int cycle = 0; cycle < settings.maxCycles; cycle++) {
    // Ouput and init
    if (printTitle)
      printSubSectionTitle((std::string) "Freeze-and-Thaw Cycle: " + (cycle + 1));
    std::vector<bool> convergedSubsystems(_activeSystems.size(), false);

    // Active system cycles
    for (unsigned int nSystem = 0; nSystem < _activeSystems.size(); nSystem++) {
      // Output and init
      auto activeSystem = _activeSystems[nSystem];

      // Keep densities on grid alive
      if (nSystem > 0 and settings.gridCutOff < 0.0) {
        for (unsigned int nSystem = 0; nSystem < _activeSystems.size(); nSystem++) {
          densOnGridControllers[nSystem] = DensityOnGridFactory<SCFMode>::produce(
              _activeSystems[nSystem]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
              supersystemgrid, 1, activeSystem->getSettings());
        }
        for (unsigned int nSystem = 0; nSystem < _passiveSystems.size(); nSystem++) {
          densOnGridControllers[nSystem + _activeSystems.size()] = DensityOnGridFactory<SCFMode>::produce(
              _passiveSystems[nSystem]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
              supersystemgrid, 1, activeSystem->getSettings());
        }
      }

      if (printTitle)
        printBigCaption((std::string) "Active System: " + activeSystem->getSystemName());
      // Store old energy of the active System
      DensityMatrix<RESTRICTED> rhoOld(activeSystem->getBasisController());
      if (activeSystem->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        rhoOld = activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix().total();
      }
      else if (activeSystem->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        rhoOld =
            activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
      }
      else {
        assert(false);
      }
      // Set up systems
      auto passiveSystems = _activeSystems;
      passiveSystems.erase(std::remove(passiveSystems.begin(), passiveSystems.end(), activeSystem), passiveSystems.end());
      passiveSystems.insert(passiveSystems.end(), _passiveSystems.begin(), _passiveSystems.end());
      // Update the active System
      if (activeSystem->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        FDETask<Options::SCF_MODES::RESTRICTED> task(activeSystem, passiveSystems);
        task.settings.embedding = settings.embedding;
        task.settings.lcSettings = settings.lcSettings;
        task.settings.gridCutOff = settings.gridCutOff;
        task.settings.smallSupersystemGrid = settings.smallSupersystemGrid;
        task.settings.finalGrid = false;
        task.settings.calculateMP2Energy = false;
        task.settings.calculateSolvationEnergy = settings.calculateSolvationEnergy;
        task.settings.firstPassiveSystemIndex = _activeSystems.size() - 1;
        task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::NONE;
        task.settings.mp2Type = settings.mp2Type;
        task.settings.skipSCF = settings.onlyFinalEnergyEvaluation;
        if (settings.gridCutOff < 0)
          task.setSuperSystemGrid(supersystemgrid);
        if (molecularSurfaceIsInitialized)
          task.settings.initializeSuperMolecularSurface = false;
        if (settings.calculateUnrelaxedMP2Density.size() == _activeSystems.size()) {
          task.settings.calculateUnrelaxedMP2Density = settings.calculateUnrelaxedMP2Density[nSystem];
        }
        task.run();
        supersystemgrid = task.getSuperSystemGrid();
      }
      else if (activeSystem->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        FDETask<Options::SCF_MODES::UNRESTRICTED> task(activeSystem, passiveSystems);
        task.settings.embedding = settings.embedding;
        task.settings.gridCutOff = settings.gridCutOff;
        task.settings.smallSupersystemGrid = settings.smallSupersystemGrid;
        task.settings.finalGrid = false;
        task.settings.calculateMP2Energy = false;
        task.settings.calculateSolvationEnergy = settings.calculateSolvationEnergy;
        task.settings.firstPassiveSystemIndex = _activeSystems.size() - 1;
        task.settings.embedding.dispersion = Options::DFT_DISPERSION_CORRECTIONS::NONE;
        task.settings.skipSCF = settings.onlyFinalEnergyEvaluation;
        if (settings.gridCutOff < 0)
          task.setSuperSystemGrid(supersystemgrid);
        if (molecularSurfaceIsInitialized)
          task.settings.initializeSuperMolecularSurface = false;
        task.run();
        supersystemgrid = task.getSuperSystemGrid();
      }
      else {
        assert(false);
      }
      molecularSurfaceIsInitialized = true;
      // Converged subsystem?
      double deltaRho = std::numeric_limits<double>::infinity();
      if (activeSystem->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        deltaRho =
            (rhoOld -
             activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix().total())
                .array()
                .square()
                .sum();
        deltaRho = sqrt(deltaRho);
      }
      else if (activeSystem->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        deltaRho =
            (rhoOld -
             activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total())
                .array()
                .square()
                .sum();
        deltaRho = sqrt(deltaRho);
      }
      else {
        assert(false);
      }
      if (fabs(deltaRho) < settings.convThresh) {
        convergedSubsystems[nSystem] = true;
      }

      // Calculate total energies
      assert((activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT ||
              activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) &&
             "Unknown electronic theory!");
      auto energyContribution = (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
                                    ? ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT
                                    : ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT;
      if (settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::ALMO)
        energyContribution = (activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT)
                                 ? ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY
                                 : ENERGY_CONTRIBUTIONS::HF_ENERGY;
      if (activeSystem->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        superSysE =
            activeSystem->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energyContribution);
      }
      else if (activeSystem->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        superSysE =
            activeSystem->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energyContribution);
      }
      else {
        assert(false);
      }
      // Output
      if (printTitle)
        printf("\n");
      if (printTitle)
        printSmallCaption("Freeze-and-Thaw Info");
      if (printTitle)
        printf("Freeze-and-Thaw Cycle               %16d\n", cycle + 1);
      if (printTitle)
        printf("Active subsystem:                   %16d\n", nSystem + 1);
      if (printTitle)
        printf("Density change (abs. deviation):        %e\n", deltaRho);
      if (printTitle)
        printf("Total supersystem energy:         %18.10f\n", superSysE);
      if (printTitle)
        printf("\n");
    }

    // Convergence check.
    if (std::find(convergedSubsystems.begin(), convergedSubsystems.end(), false) == convergedSubsystems.end()) {
      if (settings.printResults) {
        printSubSectionTitle((std::string) "Converged after " + (cycle + 1) + " freeze-and-thaw cycles!");
      }
      break;
    }
    // Accelerate convergence
    if (settings.useConvAcceleration)
      fatConvAcc->accelerateConvergence();
  }

  this->finalEnergyEvaluation();

  // allow for freeing of Libint engines
  libint.freeEngines(LIBINT_OPERATOR::kinetic, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::nuclear, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::overlap, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);

  // reset print options
  iOOptions.printSCFCycleInfo = info;
  iOOptions.printSCFResults = results;
  iOOptions.gridAccuracyCheck = gridAcc;

  // Clean up cache.
  cleanUp();
}

template<Options::SCF_MODES SCFMode>
void FreezeAndThawTask<SCFMode>::finalEnergyEvaluation() {
  calculateNonAdditiveDispersionCorrection();
  if ((settings.gridCutOff >= 0.0 and settings.finalEnergyEvaluation) or settings.onlyFinalEnergyEvaluation) {
    printSmallCaption("Creating final grid for energy evaluation");
    // atoms of all subsystems
    auto finalGridGeometry = std::make_shared<Geometry>();

    for (auto sys : _activeSystems) {
      *finalGridGeometry += *sys->getGeometry();
    }

    for (auto sys : _passiveSystems) {
      *finalGridGeometry += *sys->getGeometry();
    }
    finalGridGeometry->deleteIdenticalAtoms();

    Options::GRID_PURPOSES finalGridacc =
        (settings.smallSupersystemGrid) ? Options::GRID_PURPOSES::SMALL : Options::GRID_PURPOSES::DEFAULT;
    auto finalGrid =
        AtomCenteredGridControllerFactory::produce(finalGridGeometry, _activeSystems[0]->getSettings().grid, finalGridacc);

    auto xcfunc = settings.embedding.customNaddXCFunc.basicFunctionals.size()
                      ? Functional(settings.embedding.customNaddXCFunc)
                      : resolveFunctional(settings.embedding.naddXCFunc);
    auto kinefunc = settings.embedding.customNaddKinFunc.basicFunctionals.size()
                        ? Functional(settings.embedding.customNaddKinFunc)
                        : resolveFunctional(settings.embedding.naddKinFunc);
    bool gga(kinefunc.getFunctionalClass() != CompositeFunctionals::CLASSES::LDA or
             xcfunc.getFunctionalClass() != CompositeFunctionals::CLASSES::LDA);

    std::unique_ptr<DensityOnGrid<SCFMode>> superSysDens(new DensityOnGrid<SCFMode>(finalGrid));
    std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> superSysdensGrad(nullptr);
    if (gga)
      superSysdensGrad = std::move(makeGradientPtr<DensityOnGrid<SCFMode>>(finalGrid));

    double naddXcEnergy = 0.0;
    double naddKinEnergy = 0.0;

    for (unsigned int nSystem = 0; nSystem < _activeSystems.size(); nSystem++) {
      auto calc = DensityOnGridFactory<SCFMode>::produce(
          _activeSystems[nSystem]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(), finalGrid,
          1, _activeSystems[nSystem]->getSettings());
      *superSysDens += calc->getDensityOnGrid();
      if (gga) {
        superSysdensGrad->x += calc->getDensityGradientOnGrid().x;
        superSysdensGrad->y += calc->getDensityGradientOnGrid().y;
        superSysdensGrad->z += calc->getDensityGradientOnGrid().z;
      }

      FunctionalLibrary<SCFMode> flib(128);
      naddXcEnergy -= flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, xcfunc, calc, 0).energy;
      naddKinEnergy -= flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, kinefunc, calc, 0).energy;
    }

    for (unsigned int nSystem = 0; nSystem < _passiveSystems.size(); nSystem++) {
      auto calc = DensityOnGridFactory<SCFMode>::produce(
          _passiveSystems[nSystem]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(), finalGrid,
          1, _passiveSystems[nSystem]->getSettings());
      *superSysDens += calc->getDensityOnGrid();
      if (gga) {
        superSysdensGrad->x += calc->getDensityGradientOnGrid().x;
        superSysdensGrad->y += calc->getDensityGradientOnGrid().y;
        superSysdensGrad->z += calc->getDensityGradientOnGrid().z;
      }
      FunctionalLibrary<SCFMode> flib(128);
      naddXcEnergy -= flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, xcfunc, calc, 0).energy;
      naddKinEnergy -= flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, kinefunc, calc, 0).energy;
    }

    auto superSysDensOnGrid =
        std::make_shared<CustomDensityOnGridController<SCFMode>>(std::move(superSysDens), std::move(superSysdensGrad));
    FunctionalLibrary<SCFMode> flib(128);
    naddXcEnergy += flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, xcfunc, superSysDensOnGrid, 0).energy;
    naddKinEnergy += flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS, kinefunc, superSysDensOnGrid, 0).energy;

    for (auto& sys : _activeSystems) {
      auto eCont = sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
      eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_XC, naddXcEnergy);
      eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_KINETIC, naddKinEnergy);
    }
  }
  OutputControl::nOut << std::endl;
  if (settings.printResults) {
    printSubSectionTitle("Final Freeze-and-Thaw Energies");
    for (unsigned int nSystem = 0; nSystem < _activeSystems.size(); nSystem++) {
      printBigCaption((std::string) "Active System: " + _activeSystems[nSystem]->getSystemName());
      auto eCont = _activeSystems[nSystem]->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
      eCont->printAllComponents();
    }
  }
}

template<Options::SCF_MODES SCFMode>
void FreezeAndThawTask<SCFMode>::cleanUp() {
  for (auto sys : _activeSystems) {
    if (not settings.keepCoulombCache) {
      std::string fileName = sys->getSettings().path + sys->getSystemName() + ".pasCoulomb.h5";
      std::remove(fileName.c_str());
    }
    // TODO: Remove this as soon as the destruction of the system controllers works properly.
    if (settings.embedding.pcm.use)
      sys->template getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::FDE)->cleanUpDisk();
  }
  // TODO: Remove this as soon as the destruction of the system controllers works properly.
  if (settings.embedding.pcm.use) {
    for (auto sys : _passiveSystems) {
      sys->template getElectrostaticPotentialOnMolecularSurfaceController<SCFMode>(MOLECULAR_SURFACE_TYPES::FDE)->cleanUpDisk();
    }
  }
}

template<Options::SCF_MODES SCFMode>
void FreezeAndThawTask<SCFMode>::calculateNonAdditiveDispersionCorrection() {
  double nadDispCorrection = 0.0;
  double solvScaledNadDispCorrection = 0.0;
  auto activeSystem = _activeSystems[0];
  if (settings.embedding.dispersion != Options::DFT_DISPERSION_CORRECTIONS::NONE) {
    Timings::takeTime("FDE -    Non-Add. Disper.");
    auto environmentGeometry = std::make_shared<Geometry>();
    for (unsigned int iAct = 1; iAct < _activeSystems.size(); ++iAct)
      *environmentGeometry += *_activeSystems[iAct]->getGeometry();
    for (auto sys : _passiveSystems)
      *environmentGeometry += *sys->getGeometry();
    environmentGeometry->deleteIdenticalAtoms();
    if (settings.calculateSolvationEnergy) {
      double scaling = 0.75;
      solvScaledNadDispCorrection = scaling * DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
                                                  settings.embedding.dispersion, activeSystem->getGeometry(),
                                                  environmentGeometry, settings.embedding.naddXCFunc);
    }
    if (settings.finalEnergyEvaluation) {
      auto supersystemGeometry = std::make_shared<Geometry>();
      *supersystemGeometry += *activeSystem->getGeometry();
      *supersystemGeometry += *environmentGeometry;
      supersystemGeometry->deleteIdenticalAtoms();
      nadDispCorrection += DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
          settings.embedding.dispersion, supersystemGeometry, settings.embedding.naddXCFunc);
      for (auto sys : _passiveSystems)
        nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
            settings.embedding.dispersion, sys->getGeometry(), settings.embedding.naddXCFunc);
      for (auto sys : _activeSystems)
        nadDispCorrection -= DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
            settings.embedding.dispersion, sys->getGeometry(), settings.embedding.naddXCFunc);
    } // if settings.finalEnergyEvaluation
    Timings::timeTaken("FDE -    Non-Add. Disper.");
  }
  if (settings.calculateSolvationEnergy) {
    auto eCont = activeSystem->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_NAD_DISP, solvScaledNadDispCorrection);
  }
  for (auto& sys : _activeSystems) {
    auto eCont = sys->template getElectronicStructure<SCFMode>()->getEnergyComponentController();
    eCont->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::FDE_NAD_DISP, nadDispCorrection);
  }
}

template class FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>;
template class FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
