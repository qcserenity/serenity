/**
 * @file   NumericalGeomGradCalc.cpp
 *
 * @date   May 18, 2015
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
#include "geometry/gradients/NumericalGeomGradCalc.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "energies/EnergyContributions.h"
#include "geometry/Atom.h"
#include "geometry/AtomType.h"
#include "geometry/Geometry.h"
#include "math/Matrix.h"
#include "parameters/Constants.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/ScfTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
NumericalGeomGradCalc<SCFMode>::NumericalGeomGradCalc(double stepsize) : _delta(stepsize) {
}

template<Options::SCF_MODES SCFMode>
void NumericalGeomGradCalc<SCFMode>::calcFDEGradients(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                      std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                      CompositeFunctionals::KINFUNCTIONALS naddKinFunc,
                                                      CompositeFunctionals::XCFUNCTIONALS naddXCFunc,
                                                      double FDEgridCutOff, int FaTmaxCycles, double FaTenergyConvThresh,
                                                      Options::DFT_DISPERSION_CORRECTIONS dispersion) {
  iOOptions.printGridInfo = false;
  const auto oldPrintLevel = GLOBAL_PRINT_LEVEL;
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  FreezeAndThawTask<SCFMode> task(activeSystems, environmentSystems);
  task.settings.embedding.naddKinFunc = naddKinFunc;
  task.settings.embedding.naddXCFunc = naddXCFunc;
  task.settings.gridCutOff = FDEgridCutOff;
  task.settings.maxCycles = FaTmaxCycles;
  task.settings.convThresh = FaTenergyConvThresh;
  task.settings.embedding.dispersion = dispersion;
  task.run();
  GLOBAL_PRINT_LEVEL = oldPrintLevel;

  for (unsigned int i = 0; i < activeSystems.size(); i++) {
    Matrix<double> gradients(activeSystems[i]->getGeometry()->getNAtoms(), 3);
    gradients.setZero();
    auto electronicStructure = activeSystems[i]->getElectronicStructure<SCFMode>();

    unsigned int count = 0;
    for (auto atom : activeSystems[i]->getGeometry()->getAtoms()) {
      ++count;
      /*
       * X
       */
      atom->addToX(_delta);
      printf("%4s %3d %2s: X + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 0) = electronicStructure->getEnergy();
      atom->addToX(-2.0 * _delta);
      printf("%4s %3d %2s: X - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 0) -= electronicStructure->getEnergy();

      gradients(count - 1, 0) /= (2.0 * _delta);
      atom->addToX(_delta);
      /*
       * Y
       */
      atom->addToY(_delta);
      printf("%4s %3d %2s: Y + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 1) = electronicStructure->getEnergy();

      atom->addToY(-2.0 * _delta);
      printf("%4s %3d %2s: Y - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 1) -= electronicStructure->getEnergy();

      gradients(count - 1, 1) /= (2.0 * _delta);
      atom->addToY(_delta);
      /*
       * Z
       */
      atom->addToZ(_delta);
      printf("%4s %3d %2s: Z + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 2) = electronicStructure->getEnergy();

      atom->addToZ(-2.0 * _delta);
      printf("%4s %3d %2s: Z - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
      task.run();
      gradients(count - 1, 2) -= electronicStructure->getEnergy();

      gradients(count - 1, 2) /= (2.0 * _delta);
      atom->addToZ(_delta);
      /*
       *  Store
       */
    }

    activeSystems[i]->getGeometry()->setGradients(gradients);
    /*
     * Reset
     * TODO what is the use of this? Why calculate the initial electronic structure again?
     * Has the active electronic structure even been changed in the system?
     */
    task.run();
  }
}

template<Options::SCF_MODES SCFMode>
void NumericalGeomGradCalc<SCFMode>::calcGradients(std::shared_ptr<SystemController> systemController) {
  iOOptions.printGridInfo = false;
  ScfTask<SCFMode> scf(systemController);
  scf.settings.restart = false;
  auto electronicStructure = systemController->getElectronicStructure<SCFMode>();
  Matrix<double> gradients(systemController->getGeometry()->getNAtoms(), 3);
  gradients.setZero();

  unsigned int count = 0;
  for (auto atom : systemController->getGeometry()->getAtoms()) {
    ++count;
    /*
     * X
     */
    atom->addToX(_delta);
    printf("%4s %3d %2s: X + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 0) = electronicStructure->getEnergy();
    atom->addToX(-2.0 * _delta);
    printf("%4s %3d %2s: X - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 0) -= electronicStructure->getEnergy();

    gradients(count - 1, 0) /= (2.0 * _delta);
    atom->addToX(_delta);
    /*
     * Y
     */
    atom->addToY(_delta);
    printf("%4s %3d %2s: Y + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 1) = electronicStructure->getEnergy();

    atom->addToY(-2.0 * _delta);
    printf("%4s %3d %2s: Y - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 1) -= electronicStructure->getEnergy();

    gradients(count - 1, 1) /= (2.0 * _delta);
    atom->addToY(_delta);
    /*
     * Z
     */
    atom->addToZ(_delta);
    printf("%4s %3d %2s: Z + %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 2) = electronicStructure->getEnergy();

    atom->addToZ(-2.0 * _delta);
    printf("%4s %3d %2s: Z - %7.5f \n", "", count, atom->getAtomType()->getElementSymbol().c_str(), _delta * BOHR_TO_ANGSTROM);
    scf.run();
    gradients(count - 1, 2) -= electronicStructure->getEnergy();

    gradients(count - 1, 2) /= (2.0 * _delta);
    atom->addToZ(_delta);
    /*
     *  Store
     */
  }

  systemController->getGeometry()->setGradients(gradients);
  /*
   * Reset
   * TODO what is the use of this? Why calculate the initial electronic structure again?
   * Has the active electronic structure even been changed in the system?
   */
  scf.run();
}

template class NumericalGeomGradCalc<Options::SCF_MODES::RESTRICTED>;
template class NumericalGeomGradCalc<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
