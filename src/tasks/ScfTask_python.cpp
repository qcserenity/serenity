/**
 * @file ScfTask_python.cpp
 *
 * @date Apr 25, 2016
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
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;
void export_ScfTask(py::module& spy) {
  py::class_<ScfTaskSettings>(spy, "ScfTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("restart", &ScfTaskSettings::restart)
      .def_readwrite("mp2Type", &ScfTaskSettings::mp2Type)
      .def_readwrite("maxResidual", &ScfTaskSettings::maxResidual)
      .def_readwrite("maxCycles", &ScfTaskSettings::maxCycles)
      .def_readwrite("skipSCF", &ScfTaskSettings::skipSCF)
      .def_readwrite("allowNotConverged", &ScfTaskSettings::allowNotConverged)
      .def_readwrite("calculateMP2Energy", &ScfTaskSettings::calculateMP2Energy)
      .def_readwrite("exca", &ScfTaskSettings::exca)
      .def_readwrite("excb", &ScfTaskSettings::excb)
      .def_readwrite("momCycles", &ScfTaskSettings::momCycles)
      .def_readwrite("lcSettings", &ScfTaskSettings::lcSettings);

  py::class_<ScfTask<Options::SCF_MODES::RESTRICTED>>(spy, "ScfTask_R")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &ScfTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &ScfTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &ScfTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<ScfTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "ScfTask_U")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &ScfTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &ScfTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &ScfTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}