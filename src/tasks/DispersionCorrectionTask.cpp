/**
 * @file   DispersionCorrectionTask.cpp
 *
 * @date   Dec 02, 2015
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
#include "tasks/DispersionCorrectionTask.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "io/FormattedOutput.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

DispersionCorrectionTask::DispersionCorrectionTask(std::shared_ptr<SystemController> systemController)
  : _systemController(systemController) {
}

void DispersionCorrectionTask::run() {
  printSubSectionTitle("DFT Dispersion Correction Calculation: " + _systemController->getSystemName());

  printSmallCaption("Energy Correction");
  /*
   * If no XC-functional was specified in the DispersionCorrectionTask, use
   * the XC-functional specified in the system.
   */
  if (settings.functional == CompositeFunctionals::XCFUNCTIONALS::NONE) {
    settings.functional = (_systemController->getSettings()).dft.functional;
  }
  auto energy = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      settings.dispType, _systemController->getGeometry(), settings.functional);
  /*
   * Also print for which XC-functional the correction is applied
   */
  std::string functional;
  auto func = settings.functional;
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(functional, func);
  printf("%4s Energy Correction (%15s): %+15.10f Hartree \n", "", functional.c_str(), energy);
  printf("\n \n");
  if (settings.gradient) {
    printSmallCaption("Gradient Correction");
    auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection(
        settings.dispType, _systemController->getGeometry(), settings.functional);
    printf("%4s Gradient Correction (in a.u.)\n", "");
    for (unsigned int i = 0; i < grad.rows(); ++i) {
      printf("%4s  %5d %2s %+15.10f  %+15.10f  %+15.10f\n", "", (i + 1),
             _systemController->getGeometry()->getAtoms()[i]->getAtomType()->getElementSymbol().c_str(), grad(i, 0),
             grad(i, 1), grad(i, 2));
    }
    printf("\n \n");
  }
  if (settings.hessian) {
    printSmallCaption("Hessian Correction");
    auto hess = DispersionCorrectionCalculator::calcDispersionHessianCorrection(
        settings.dispType, _systemController->getGeometry(), settings.functional);

    printf("%4s Hesian Correction (in a.u.) \n      (order xx,xy,xz,yy,yz,zz)\n \n", "");

    for (auto dir : hess) {
      for (unsigned int i = 0; i < dir.rows(); ++i) {
        printf("%4s  %5d %2s %+15.10f", "", (i + 1),
               _systemController->getGeometry()->getAtoms()[i]->getAtomType()->getElementSymbol().c_str(), dir(i, 0));
        for (unsigned int j = 1; j < dir.rows(); ++j) {
          printf("  %+15.10f", dir(i, j));
          if (j == 4) {
            printf("  ...\n");
            break;
          }
          else if (j == dir.rows() - 1)
            printf("\n");
        }
      }
      printf("\n");
    }
    printf("\n \n");
  }
}

} /* namespace Serenity */
