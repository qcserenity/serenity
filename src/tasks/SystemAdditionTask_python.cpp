/**
 * @file SystemAdditionTask_python.cpp
 *
 * @date Nov 6, 2024
 * @author Lukas Lampe
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
#include "tasks/SystemAdditionTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_SystemAdditionTask(py::module& spy) {
  py::class_<SystemAdditionTaskSettings>(spy, "SystemAdditionTaskSettings",
                                         "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("checkSuperGeom", &SystemAdditionTaskSettings::checkSuperGeom)
      .def_readwrite("checkSuperCharge", &SystemAdditionTaskSettings::checkSuperCharge)
      .def_readwrite("addOccupiedOrbitals", &SystemAdditionTaskSettings::addOccupiedOrbitals);

  py::class_<SystemAdditionTask<Options::SCF_MODES::RESTRICTED>>(spy, "SystemAdditionTask_R")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &SystemAdditionTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &SystemAdditionTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &SystemAdditionTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "SystemAdditionTask_U")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}