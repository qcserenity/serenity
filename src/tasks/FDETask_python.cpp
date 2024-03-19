/**
 * @file FDETask_python.cpp
 *
 * @date Apr 28, 2016
 * @author: Jan Unsleber
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
#include "system/SystemController.h"
#include "tasks/FDETask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_FDETask(py::module& spy) {
  py::class_<FDETaskSettings>(spy, "FDETaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("gridCutOff", &FDETaskSettings::gridCutOff)
      .def_readwrite("smallSupersystemGrid", &FDETaskSettings::smallSupersystemGrid)
      .def_readwrite("finalGrid", &FDETaskSettings::finalGrid)
      .def_readwrite("calculateUnrelaxedMP2Density", &FDETaskSettings::calculateUnrelaxedMP2Density)
      .def_readwrite("calculateMP2Energy", &FDETaskSettings::calculateMP2Energy)
      .def_readwrite("maxResidual", &FDETaskSettings::maxResidual)
      .def_readwrite("maxCycles", &FDETaskSettings::maxCycles)
      .def_readwrite("calculateEnvironmentEnergy", &FDETaskSettings::calculateEnvironmentEnergy)
      .def_readwrite("mp2Type", &FDETaskSettings::mp2Type)
      .def_readwrite("calculateSolvationEnergy", &FDETaskSettings::calculateSolvationEnergy)
      .def_readwrite("skipSCF", &FDETaskSettings::skipSCF)
      .def_readwrite("embedding", &FDETaskSettings::embedding)
      .def_readwrite("lcSettings", &FDETaskSettings::lcSettings)
      .def_readwrite("loc", &FDETaskSettings::loc);

  py::class_<FDETask<Options::SCF_MODES::RESTRICTED>>(spy, "FDETask_R")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &FDETask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &FDETask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &FDETask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<FDETask<Options::SCF_MODES::UNRESTRICTED>>(spy, "FDETask_U")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &FDETask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &FDETask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &FDETask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}