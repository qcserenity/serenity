/**
 *  @file ReadOrbitalsTask_python.cpp
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
#include "tasks/ReadOrbitalsTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

void export_ReadOrbitalsTask(py::module& spy) {
  py::class_<ReadOrbitalsTaskSettings>(spy, "ReadOrbitalsTaskSettings",
                                       "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("fileFormat", &ReadOrbitalsTaskSettings::fileFormat)
      .def_readwrite("path", &ReadOrbitalsTaskSettings::path)
      .def_readwrite("resetCoreOrbitals", &ReadOrbitalsTaskSettings::resetCoreOrbitals)
      .def_readwrite("replaceInFile", &ReadOrbitalsTaskSettings::replaceInFile);

  py::class_<ReadOrbitalsTask<Options::SCF_MODES::RESTRICTED>>(spy, "ReadOrbitalsTask_R")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &ReadOrbitalsTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &ReadOrbitalsTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &ReadOrbitalsTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<ReadOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "ReadOrbitalsTask_U")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &ReadOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &ReadOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &ReadOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
