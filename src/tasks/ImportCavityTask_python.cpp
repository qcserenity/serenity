/**
 * @file ImportCavityTask_python.cpp
 *
 * @date Dec 02, 2024
 * @author: Lukas Paetow
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
#include "tasks/ImportCavityTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;
void export_ImportCavityTask(py::module& spy) {
  py::class_<ImportCavityTaskSettings>(spy, "ImportCavityTaskSettings",
                                       "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("fdecavity", &ImportCavityTaskSettings::fdecavity)
      .def_readwrite("cavityPath", &ImportCavityTaskSettings::cavityPath)
      .def_readwrite("vdwcavityPath", &ImportCavityTaskSettings::vdwcavityPath);

  py::class_<ImportCavityTask>(spy, "ImportCavityTask")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &ImportCavityTask::run)
      .def_readwrite("settings", &ImportCavityTask::settings)
      .def_readwrite("generalSettings", &ImportCavityTask::generalSettings);
}