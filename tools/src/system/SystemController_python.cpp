/**
 * @file SystemController_python.cpp
 *
 * @date Apr 25, 2016
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

/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "geometry/Geometry.h"
#include "grid/GridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "scf/initialGuess/InitialGuessFactory.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

std::shared_ptr<SystemController> createSystemControllerPtr1(Settings& settings) {
  return std::make_shared<SystemController>(settings);
}

std::shared_ptr<SystemController> add(std::shared_ptr<SystemController> sys1, std::shared_ptr<SystemController> sys2) {
  return (*sys1) + (*sys2);
}

void guess(std::shared_ptr<SystemController> system) {
  if (!system->hasElectronicStructure<RESTRICTED>() and system->getSettings().scfMode == RESTRICTED) {
    system->setElectronicStructure(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>(
        InitialGuessFactory::produce<Options::SCF_MODES::RESTRICTED>(system->getSettings().scf.initialguess)
            ->calculateInitialGuess(system)));
  }
  if (!system->hasElectronicStructure<UNRESTRICTED>() and system->getSettings().scfMode == UNRESTRICTED) {
    system->setElectronicStructure(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>(
        InitialGuessFactory::produce<Options::SCF_MODES::UNRESTRICTED>(system->getSettings().scf.initialguess)
            ->calculateInitialGuess(system)));
  }
}

double getEnergy1(std::shared_ptr<SystemController> system, Options::SCF_MODES mode, ENERGY_CONTRIBUTIONS energy) {
  if (mode == Options::SCF_MODES::RESTRICTED) {
    return system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energy);
  }
  else {
    return system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energy);
  }
}

double getEnergy2(std::shared_ptr<SystemController> system, ENERGY_CONTRIBUTIONS energy) {
  if (system->getLastSCFMode() == Options::SCF_MODES::RESTRICTED) {
    return system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(energy);
  }
  else {
    return system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy(energy);
  }
}

double getEnergy3(std::shared_ptr<SystemController> system, Options::SCF_MODES mode) {
  if (mode == Options::SCF_MODES::RESTRICTED) {
    return system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  }
  else {
    return system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy();
  }
}

double getEnergy4(std::shared_ptr<SystemController> system) {
  if (system->getLastSCFMode() == Options::SCF_MODES::RESTRICTED) {
    return system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  }
  else {
    return system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergy();
  }
}

int getNEl(std::shared_ptr<SystemController> system) {
  return system->getNElectrons<RESTRICTED>();
}

Eigen::MatrixXd overlap(std::shared_ptr<SystemController> system) {
  return system->getOneElectronIntegralController()->getOverlapIntegrals();
}

std::shared_ptr<SystemController> getSptr(std::shared_ptr<SystemController> system) {
  return system->getSharedPtr();
}

void export_SystemController(py::module& spy) {
  spy.def("combine", &add);
  py::class_<SystemController, std::shared_ptr<SystemController>>(spy, "System")
      .def(py::init(&createSystemControllerPtr1), "@brief The basic constructor for a system\n"
                                                  "\n"
                                                  "@param libserenipy.Settings The Settings for this System.")
      .def(py::init<std::shared_ptr<Geometry>, Settings>())
      .def("get", &getSptr)
      .def("guess", &guess)
      .def("getEnergy", &getEnergy1,
           "@brief Returns a specific energy for the last run SCF of the given mode.\n"
           "\n"
           "@param SCF_MODES The mode [RESTRICTED, UNRESTRICTED].\n"
           "@param ENERGY_CONTRIBUTIONS The energy contribution to be returned.\n"
           "\n"
           "@returns The energy.")
      .def("getEnergy", &getEnergy2,
           "@brief Returns a specific energy for the last run SCF.\n"
           "\n"
           "@param SCF_MODE The mode [RESTRICTED, UNRESTRICTED].\n"
           "\n"
           "@returns The energy.")
      .def("getEnergy", &getEnergy3,
           "@brief Returns the total energy for the last run SCF of the given mode.\n"
           "\n"
           "@param ENERGY_CONTRIBUTIONS The energy contribution to be returned.\n"
           "\n"
           "@returns The energy.")
      .def("getEnergy", &getEnergy4,
           "@brief Returns the total energy of the last run SCF.\n"
           "\n"
           "@returns The energy.")
      .def("getGeometry", &SystemController::getGeometry,
           "@brief Returns the current geometry.\n"
           "\n"
           "@returns The geometry.")
      .def("getOverlap", &overlap)
      .def("getSpin", &SystemController::getSpin)
      .def("getCharge", &SystemController::getCharge)
      .def("getNEletrons", &getNEl)
      .def("getBasis", &SystemController::getBasisController, py::arg("basisPurpose") = Options::BASIS_PURPOSES::DEFAULT)
      .def("getGrid", &SystemController::getGridController, py::arg("gridPurpose") = Options::GRID_PURPOSES::DEFAULT)
      .def("getElectronicStructure_R", &SystemController::getElectronicStructure<Options::SCF_MODES::RESTRICTED>,
           "@brief Returns the current restricted electronic structure.\n"
           "       If there is none, a SCF will be run to generate one.\n"
           "@returns The restricted electronic structure.")
      .def("getElectronicStructure_U", &SystemController::getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>,
           "@brief Returns the current unrestricted electronic structure.\n"
           "       If there is none, a SCF will be run to generate one.\n"
           "@returns The unrestricted electronic structure.")
      .def("getSystemName", &SystemController::getSystemName,
           "@returns the name of the controlled molecular system. It should be unique.")
      .def("getSettings", &SystemController::getSettings,
           "@returns the underlying configuration. All data received from this SystemController is\n"
           "constructed by using this configuration.");
}
