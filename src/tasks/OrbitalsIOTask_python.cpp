/**
 *  @file OrbitalsIOTask_python.cpp
 *
 *  @date   Sep 13, 2022
 *  @author Moritz Bensberg
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
#include "tasks/OrbitalsIOTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

void export_OrbitalsIOTask(py::module& spy) {
  py::class_<OrbitalsIOTaskSettings>(spy, "OrbitalsIOTaskSettings",
                                     "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("fileFormat", &OrbitalsIOTaskSettings::fileFormat)
      .def_readwrite("path", &OrbitalsIOTaskSettings::path)
      .def_readwrite("resetCoreOrbitals", &OrbitalsIOTaskSettings::resetCoreOrbitals)
      .def_readwrite("replaceInFile", &OrbitalsIOTaskSettings::replaceInFile)
      .def_readwrite("write", &OrbitalsIOTaskSettings::write);

  py::class_<OrbitalsIOTask<Options::SCF_MODES::RESTRICTED>>(spy, "OrbitalsIOTask_R")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &OrbitalsIOTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &OrbitalsIOTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &OrbitalsIOTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<OrbitalsIOTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "OrbitalsIOTask_U")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &OrbitalsIOTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &OrbitalsIOTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &OrbitalsIOTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
