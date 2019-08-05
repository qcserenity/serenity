/**
 * @file   ScfTask.cpp
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#include "tasks/ScfTask.h"
/* Include Serenity Internal Headers */
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "postHF/MPn/MP2.h"
#include "data/OrbitalController.h"
#include "potentials/bundles/PotentialBundle.h"
#include "postHF/MPn/RIMP2.h"
#include "scf/SCFAnalysis.h"
#include "scf/Scf.h"
#include "misc/SerenityError.h"
#include "system/SystemController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {
using namespace std;

template<Options::SCF_MODES SCFMode>ScfTask<SCFMode>::ScfTask(
      const shared_ptr<SystemController> systemController) :
    _systemController(systemController) {
  assert(_systemController);
}


template<Options::SCF_MODES SCFMode>void ScfTask<SCFMode>::run() {

  if (this->settings.restart) {
    std::cout << _systemController->getSettings().load << std::endl;
    if(_systemController->getSettings().load.empty()){
      throw SerenityError(
      "Option restart in SCFTask only available in combination with option load in system block!");
    }
    shared_ptr<OrbitalController<SCFMode> > orbitals;
    try{
      orbitals=std::make_shared<OrbitalController<SCFMode>>(_systemController->getSettings().load+"tmp",
          _systemController->getBasisController(),_systemController->getSettings().identifier);
    }catch(...){
      std::cout << "No temporary orbital files found. Looking for converged orbital files";
      orbitals=std::make_shared<OrbitalController<SCFMode>>(_systemController->getSettings().load+_systemController->getSettings().name,
          _systemController->getBasisController(),_systemController->getSettings().identifier);
    }
    _systemController->setElectronicStructure(
        std::make_shared<ElectronicStructure<SCFMode> >(
            orbitals,
            _systemController->getOneElectronIntegralController(),
            _systemController->getNOccupiedOrbitals<SCFMode>()));
  }

  bool guess = !(_systemController->hasElectronicStructure<SCFMode>());

  auto& settings = _systemController->getSettings();
  std::shared_ptr<PotentialBundle<SCFMode> > potentials;
  if (settings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
    potentials = _systemController->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  } else if (settings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    potentials = _systemController->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
        Options::GRID_PURPOSES::SMALL);
  } else {
    std::cout << "ERROR: None existing electronicStructureTheory requested." << std::endl;
    assert(false);
  }
  auto es = _systemController->getElectronicStructure<SCFMode>();

  if(this->settings.fractionalDegeneracy){
    es->getDensityMatrixController()->setFractionalDegeneracy(true);
  }



  // Output of the main options
  if (iOOptions.printSCFResults){
    printSubSectionTitle("Main SCF Options");
    auto m = settings.method;
    std::string method;
    Options::resolve<Options::ELECTRONIC_STRUCTURE_THEORIES>(method,m);
    printf("%4s Method:                %15s\n","",method.c_str());\
    if (settings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
      std::string functional;
      auto func = settings.dft.functional;
      Options::resolve<Options::XCFUNCTIONALS>(functional,func);
      printf("%4s Functional:            %15s\n","",functional.c_str());
      printf("%4s Grid Accuracy:         %13d/%1d\n","",settings.grid.smallGridAccuracy,settings.grid.accuracy);
      std::string fitting;
      auto fit = settings.dft.densityFitting;
      Options::resolve<Options::DENS_FITS>(fitting,fit);
      printf("%4s Density Fitting:       %15s\n","",fitting.c_str());
      std::string dispersion;
      auto disp = settings.dft.dispersion;
      Options::resolve<Options::DFT_DISPERSION_CORRECTIONS>(dispersion,disp);
      printf("%4s Dispersion Correction: %15s\n","",dispersion.c_str());
    }
    printf("%4s Basis Set:             %15s\n","",settings.basis.label.c_str());
    if(_systemController->getGeometry()->hasAtomsWithECPs()) {
      printf("%4s ECP Start:             %15d\n","",settings.basis.firstECP);
    }
    if(guess){
      auto ig = settings.scf.initialguess;
      std::string init_guess;
      Options::resolve<Options::INITIAL_GUESSES>(init_guess,ig);
      printf("%4s Initial Guess:         %15s\n","",init_guess.c_str());
    }
    printSubSectionTitle("SCF");
  }

  auto energyComponentController = es->getEnergyComponentController();
  if (settings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
    // disable orbital output for dft
    auto tmp1 = iOOptions.printFinalOrbitalEnergies;
    auto tmp2 = iOOptions.gridAccuracyCheck;
    iOOptions.printFinalOrbitalEnergies = false;
    iOOptions.gridAccuracyCheck = false;

    Scf<SCFMode>::perform(settings,es,potentials);

    iOOptions.printFinalOrbitalEnergies = tmp1;
    iOOptions.gridAccuracyCheck = tmp2;
    // run final cycle on larger (default) grid
    potentials = _systemController->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
        Options::GRID_PURPOSES::DEFAULT);

    auto orbitalController =  es->getMolecularOrbitals();
    orbitalController->updateOrbitals(potentials->getFockMatrix(es->getDensityMatrix(),energyComponentController),
                                     es->getOneElectronIntegralController());

    //check for double hybrid functional
    auto functional = FunctionalClassResolver::resolveFunctional(settings.dft.functional);
    double MP2Correlation = 0.0;
    if (functional.isDoubleHybrid()){
      //perform MP2 for double hybrids
      if (settings.dft.densityFitting == Options::DENS_FITS::RI){
        RIMP2<SCFMode> rimp2(_systemController, functional.getssScaling(), functional.getosScaling());
        MP2Correlation = rimp2.calculateCorrection();
      } else {
        assert(SCFMode==RESTRICTED && "MP2 is not available for unrestricted systems please use RI-MP2.");
        MP2EnergyCorrector<RESTRICTED> mp2EnergyCorrector(_systemController, functional.getssScaling(), functional.getosScaling());
        MP2Correlation = mp2EnergyCorrector.calculateElectronicEnergy();
      }
      MP2Correlation *= functional.getHfCorrelRatio();
    }

    // add energy to the EnergyController
    energyComponentController->addOrReplaceComponent(
        std::pair<ENERGY_CONTRIBUTIONS,double>(ENERGY_CONTRIBUTIONS::KS_DFT_PERTURBATIVE_CORRELATION,MP2Correlation));

  } else {
    Scf<SCFMode>::perform(settings,es,potentials);
  }


  if (settings.method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
      auto dispCorrection=DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
          settings.dft.dispersion,
          _systemController->getGeometry(),
          settings.dft.functional);
      es->getEnergyComponentController()->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION,dispCorrection);
  }

  /*
   * Print final results
   */
  if (iOOptions.printSCFResults){
    printSubSectionTitle("Final SCF Results");
    energyComponentController->printAllComponents();
    if (iOOptions.printFinalOrbitalEnergies){
      print("");
      printSmallCaption("Orbital Energies");
      es->printMOEnergies();
      print("");
    }
    //SCF Analysis
    printSmallCaption("Additional Analysis");
      SCFAnalysis<SCFMode> scfAn(
          _systemController,
          es->getOneElectronIntegralController(),
          energyComponentController);
      auto s2val = scfAn.S2();
      auto virialRatio = scfAn.VirialRatio();
      printf("\n   -<V>/<T> = %4.3f ", virialRatio);
      printf("\n      <S*S> = %4.3f ", s2val);
      double S=fabs(0.5 * _systemController->getSpin());
      printf("\n    S*(S+1) = %4.3f ", S*(S+1));
      printf("\n          C = %4.3f \n\n", s2val - S*(S+1));
  }

  return;
}


template class ScfTask<Options::SCF_MODES::RESTRICTED>;
template class ScfTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
