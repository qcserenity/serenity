/**
 * @file   GradientTask.cpp
 *
 * @date   Mar 23, 2015
 * @author Kevin Klahr
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
#include "tasks/GradientTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "data/ElectronicStructure.h"
#include "tasks/FDETask.h"
#include "io/FormattedOutput.h"
#include "tasks/FreezeAndThawTask.h"
#include "geometry/Geometry.h"
#include "geometry/gradients/GeometryGradientCalculator.h"
#include "io/IOOptions.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "potentials/NAddFuncPotential.h"
#include "geometry/gradients/NumericalGeomGradCalc.h"
#include "data/OrbitalController.h"
#include "system/SystemController.h"
#include "misc/Timing.h"


namespace Serenity {
using namespace std;


template<Options::SCF_MODES SCFMode>
GradientTask<SCFMode>::GradientTask(
    const std::vector<std::shared_ptr<SystemController> >& activeSystems,
    const std::vector<std::shared_ptr<SystemController> >& passiveSystems) :
    		    _activeSystems(activeSystems),
    		    _passiveSystems(passiveSystems){
  assert(_activeSystems[0]);
}

template<Options::SCF_MODES SCFMode>
GradientTask<SCFMode>::~GradientTask() {
}

template<Options::SCF_MODES SCFMode>
void GradientTask<SCFMode>::run() {

  takeTime("Gradient Calculation");
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool check (iOOptions.gridAccuracyCheck =false);
  int timings(iOOptions.timingsPrintLevel);
  iOOptions.gridAccuracyCheck =false;
  if (!(settings.print)){
    iOOptions.printSCFCycleInfo =false;
    iOOptions.printSCFResults =false;
    iOOptions.timingsPrintLevel = 0;
  }

  if (_passiveSystems.size()==0 && _activeSystems.size()==1){
    if(settings.gradType == Options::GRADIENT_TYPES::NUMERICAL){
      NumericalGeomGradCalc<SCFMode> numGradCalc(settings.numGradStepSize);
      if (settings.print)  printSubSectionTitle("Numerical Gradient Calculation");

      numGradCalc.calcGradients(_activeSystems[0]);

      if(settings.transInvar){
        _activeSystems[0]->getGeometry()->makeGradientsTranslationallyInvariant();
      }

      _activeSystems[0]->getGeometry()->printGradients();
    } else {

      auto es = _activeSystems[0]->getElectronicStructure<SCFMode>();
      if (!es->potentialsAvailable()){
        ScfTask<SCFMode> scf(_activeSystems[0]);
        scf.run();
      }
      auto potBundle = es->getPotentialBundle();
      auto potentialGradients = potBundle->getGradients();
      Matrix<double> dispCorr(_activeSystems[0]->getGeometry()->getAtoms().size(),3);
      dispCorr.setZero();
      if(!(settings.dispersion == Options::DFT_DISPERSION_CORRECTIONS::NONE)){
      //Dispersion Correction components
      dispCorr += DispersionCorrectionCalculator::calcDispersionGradientCorrection(settings.dispersion,
          _activeSystems[0]->getGeometry(),_activeSystems[0]->getSettings().dft.functional);
      }
      Matrix<double> gradient = potentialGradients + dispCorr;
      _activeSystems[0]->getGeometry()->setGradients(gradient);
      if(settings.transInvar){
        _activeSystems[0]->getGeometry()->makeGradientsTranslationallyInvariant();
      }

      /*
       * Output
       */
      if (settings.print)  _activeSystems[0]->getGeometry()->printGradients();
    }

  } else {

    /*
     * Initial FaT
     */

    FreezeAndThawTask<SCFMode> task(_activeSystems,_passiveSystems);
    task.settings.embedding.naddKinFunc = settings.naddKinFunc;
    task.settings.embedding.naddXCFunc = settings.naddXCFunc;
    task.settings.gridCutOff = settings.FDEgridCutOff;
    task.settings.maxCycles = settings.FaTmaxCycles;
    task.settings.convThresh = settings.FaTenergyConvThresh;
    task.settings.printLevel = 0;
    task.settings.embedding.dispersion = settings.dispersion;
    task.run();

    if(settings.gradType == Options::GRADIENT_TYPES::NUMERICAL){
           NumericalGeomGradCalc<SCFMode> numGradCalc(settings.numGradStepSize);
           if (settings.print)  printSubSectionTitle("Numerical FDE Gradient Calculation");

           numGradCalc.calcFDEGradients(_activeSystems,_passiveSystems,
               settings.naddKinFunc,
               settings.naddXCFunc,
               settings.FDEgridCutOff,
               settings.FaTmaxCycles,
               settings.FaTenergyConvThresh,
               settings.dispersion);

           for (unsigned int i=0; i < _activeSystems.size(); i++){

             printBigCaption((std::string)"Active System: "+(i+1));
             if(settings.transInvar){
               _activeSystems[i]->getGeometry()->makeGradientsTranslationallyInvariant();
             }
             _activeSystems[i]->getGeometry()->printGradients();
           }

         }  else {

    // Active system cycles
    for (unsigned int nSystem=0; nSystem<_activeSystems.size();nSystem++){

      // Output and init
      printBigCaption((std::string)"Active System: "+(nSystem+1));
      auto activeSystem =  _activeSystems[nSystem];

      // Set up systems
      auto passiveSystems = _activeSystems;
      passiveSystems.erase(std::remove(passiveSystems.begin(), passiveSystems.end(), activeSystem), passiveSystems.end());

      passiveSystems.insert(passiveSystems.end(), _passiveSystems.begin(), _passiveSystems.end());


      auto supersystemGeometry = std::make_shared<Geometry>();
      *supersystemGeometry += *activeSystem->getGeometry();
      for (auto sys : passiveSystems){
        *supersystemGeometry += *sys->getGeometry();
      }


        FDETask<SCFMode> task(activeSystem,passiveSystems);
        task.settings.embedding.naddKinFunc = settings.naddKinFunc;
        task.settings.embedding.naddXCFunc = settings.naddXCFunc;
        task.settings.gridCutOff = settings.FDEgridCutOff;
        task.settings.embedding.dispersion = settings.dispersion;
        task.run();

        takeTime("Gradient Step");
        auto es = activeSystem->getElectronicStructure<SCFMode>();
        auto potBundle = es->getPotentialBundle();
        auto potentialGradients = potBundle->getGradients();
        Matrix<double> ccRepDerivative(activeSystem->getGeometry()->getAtoms().size(),3);
        ccRepDerivative.setZero();
        for (auto sys : passiveSystems){
          ccRepDerivative += CoreCoreRepulsionDerivative::calculateDerivative(activeSystem->getAtoms(),sys->getAtoms());
        }

        Matrix<double> dispCorr(activeSystem->getGeometry()->getAtoms().size(),3);
        dispCorr.setZero();
        if(!(settings.dispersion == Options::DFT_DISPERSION_CORRECTIONS::NONE)){
        dispCorr += DispersionCorrectionCalculator::calcDispersionGradientCorrection(settings.dispersion,
            supersystemGeometry,settings.naddXCFunc).topRows( activeSystem->getGeometry()->getNAtoms());
        }

        //Add interaction parts into atoms
        Matrix<double> gradientSum = potentialGradients + ccRepDerivative + dispCorr;

        activeSystem->getGeometry()->setGradients(gradientSum);

        if(settings.transInvar){
          auto supersystem = std::make_shared<SystemController>(supersystemGeometry,activeSystem->getSettings());
          activeSystem->getGeometry()->makeGradientsTranslationallyInvariant();
        }
        /*
         * Output
         */
        activeSystem->getGeometry()->printGradients();
        timeTaken(1,"Gradient Step");
      }

    }

  }

  iOOptions.printSCFCycleInfo =info;
  iOOptions.printSCFResults =results;
  iOOptions.timingsPrintLevel = timings;
  iOOptions.gridAccuracyCheck = check;
  timeTaken(3, "Gradient Calculation");

}
template class GradientTask<Options::SCF_MODES::RESTRICTED>;
template class GradientTask<Options::SCF_MODES::UNRESTRICTED>;
}/* namespace Serenity */
