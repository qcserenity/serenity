/**
 * @file HessianTask.cpp
 *
 * @date Feb 7, 2017
 * @author Kevin Klahr
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
#include "tasks/HessianTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Atom.h"
#include "geometry/gradients/HessianCalculator.h"
#include "geometry/gradients/NumericalGeomGradCalc.h"
#include "geometry/gradients/NumericalHessianCalc.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "misc/Timing.h"
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HessianTask<SCFMode>::HessianTask(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                  std::vector<std::shared_ptr<SystemController>> passiveSystems)
  : _activeSystems(activeSystems), _passiveSystems(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void HessianTask<SCFMode>::run() {
  avoidMixedSCFModes(SCFMode, _activeSystems);
  /*
   * Initial SCF
   */
  // TODO add restart option here

  /*
   * Calc gradients silently
   * TODO This should be handled differently. It's never a good idea to just set global variables
   * somewhere.
   */
  takeTime("Hessian Calculation");
  bool info(iOOptions.printSCFCycleInfo);
  bool results(iOOptions.printSCFResults);
  bool check(iOOptions.gridAccuracyCheck = false);
  int timings(iOOptions.timingsPrintLevel);
  iOOptions.printSCFCycleInfo = false;
  iOOptions.printSCFResults = false;
  iOOptions.timingsPrintLevel = 0;
  iOOptions.gridAccuracyCheck = false;

  if (_passiveSystems.size() == 0 && _activeSystems.size() == 1) {
    /*
     * Choose gradient calculator
     */
    std::unique_ptr<HessianCalculator> hessCalc;
    switch (settings.hessType) {
      case Options::HESSIAN_TYPES::NUMERICAL:
        if (settings.gradType == Options::GRADIENT_TYPES::NUMERICAL) {
          hessCalc = std::unique_ptr<HessianCalculator>(
              new NumericalHessianCalc<SCFMode>(settings.numGradStepSize, settings.numHessStepSize, settings.printToFile));
          printSubSectionTitle("Numerical Hessian Calculation");
        }
        else {
          hessCalc = std::unique_ptr<HessianCalculator>(
              new NumericalHessianCalc<SCFMode>(0.0, settings.numHessStepSize, settings.printToFile));
          printSubSectionTitle("Seminumerical Hessian Calculation");
        }

        break;
      case (Options::HESSIAN_TYPES::ANALYTICAL):
        throw SerenityError("Error: Analytical Hessian not yet implemented!");
        break;
    }

    _activeSystems[0]->getElectronicStructure<SCFMode>();
    hessCalc->calcHessian(_activeSystems[0]);
  }
  else {
    /*
     * Choose gradient calculator
     */
    std::unique_ptr<HessianCalculator> hessCalc;
    switch (settings.hessType) {
      case Options::HESSIAN_TYPES::NUMERICAL:
        if (settings.gradType == Options::GRADIENT_TYPES::NUMERICAL) {
          hessCalc = std::unique_ptr<HessianCalculator>(
              new NumericalHessianCalc<SCFMode>(settings.numGradStepSize, settings.numHessStepSize, settings.printToFile));

          throw SerenityError(
              (std::string) "Fully numerical FaT Hessian not available! Please try a seminumerical calculation!");
        }
        else {
          hessCalc = std::unique_ptr<HessianCalculator>(
              new NumericalHessianCalc<SCFMode>(0.0, settings.numHessStepSize, settings.printToFile));
          printSubSectionTitle("Seminumerical Hessian Calculation");
        }

        break;
      case (Options::HESSIAN_TYPES::ANALYTICAL):
        throw SerenityError((std::string) "Analytical Hessian NYI!");
        break;
    }
    hessCalc->calcFaTHessian(_activeSystems, _passiveSystems, settings.embedding, settings.FaTmaxCycles,
                             settings.FaTenergyConvThresh, settings.FaTgridCutOff);
  }

  iOOptions.printSCFCycleInfo = info;
  iOOptions.printSCFResults = results;
  iOOptions.timingsPrintLevel = timings;
  iOOptions.gridAccuracyCheck = check;
  timeTaken(3, "Hessian Calculation");
}

template class HessianTask<Options::SCF_MODES::RESTRICTED>;
template class HessianTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
