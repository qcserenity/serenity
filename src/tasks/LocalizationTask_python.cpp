/**
 * @file LocalizationTask_python.cpp
 *
 * @date Sep 22, 2016
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
#include "tasks/LocalizationTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_LocalizationTask(py::module& spy) {
  py::class_<LocalizationTaskSettings>(spy, "LocalizationTaskSettings",
                                       "@brief Default constructor for Settings all set to their default values.")
      .def(py::init<>())
      .def_readwrite("locType", &LocalizationTaskSettings::locType)
      .def_readwrite("maxSweeps", &LocalizationTaskSettings::maxSweeps)
      .def_readwrite("alignExponent", &LocalizationTaskSettings::alignExponent)
      .def_readwrite("useKineticAlign", &LocalizationTaskSettings::useKineticAlign)
      .def_readwrite("splitValenceAndCore", &LocalizationTaskSettings::splitValenceAndCore)
      .def_readwrite("useEnergyCutOff", &LocalizationTaskSettings::useEnergyCutOff)
      .def_readwrite("energyCutOff", &LocalizationTaskSettings::energyCutOff)
      .def_readwrite("nCoreOrbitals", &LocalizationTaskSettings::nCoreOrbitals)
      .def_readwrite("localizeVirtuals", &LocalizationTaskSettings::localizeVirtuals)
      .def_readwrite("splitVirtuals", &LocalizationTaskSettings::splitVirtuals)
      .def_readwrite("virtualEnergyCutOff", &LocalizationTaskSettings::virtualEnergyCutOff)
      .def_readwrite("nRydbergOrbitals", &LocalizationTaskSettings::nRydbergOrbitals)
      .def_readwrite("replaceVirtuals", &LocalizationTaskSettings::replaceVirtuals);

  py::class_<LocalizationTask>(spy, "LocalizationTask", "A task that runs an orbital localization.")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>>())
      .def("run", &LocalizationTask::run)
      .def_readwrite("settings", &LocalizationTask::settings)
      .def_readwrite("generalSettings", &LocalizationTask::generalSettings);
}
