/**
 * @file LocalizationTask_python.cpp
 *
 * @date Sep 22, 2016
 * @author: Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "tasks/LocalizationTask.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

void export_LocalizationTask(py::module &spy){
  py::class_<LocalizationTaskSettings>(spy,"LocalizationTaskSettings",
    "@brief Default constructor for Settings all set to their default values.")
  .def_readwrite("locType",&LocalizationTaskSettings::locType)
  .def_readwrite("maxSweeps",&LocalizationTaskSettings::maxSweeps);

  py::class_<LocalizationTask >(spy,"LocalizationTask",
      "A task that exports properties on a cubic grid.")
            .def(py::init<std::shared_ptr<SystemController> >())
            .def("run", &LocalizationTask::run )
            .def_readwrite("settings",&LocalizationTask::settings);
}
