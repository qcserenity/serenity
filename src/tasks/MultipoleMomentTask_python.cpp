/**
 * @file MultipoleMomentTask_python.cpp
 *
 * @date Apr 06, 2017
 * @author: David Schnieders
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
#include "tasks/MultipoleMomentTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

void export_MultipoleMomentTask(py::module& spy) {
  py::class_<MultipoleMomentTaskSettings>(spy, "MultipoleMomentTaskSettings",
                                          "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("highestOrder", &MultipoleMomentTaskSettings::highestOrder)
      .def_readwrite("numerical", &MultipoleMomentTaskSettings::numerical)
      .def_readwrite("origin", &MultipoleMomentTaskSettings::origin)
      .def_readwrite("printTotal", &MultipoleMomentTaskSettings::printTotal)
      .def_readwrite("printFragments", &MultipoleMomentTaskSettings::printFragments);

  py::class_<MultipoleMomentTask>(spy, "MultipoleMomentTask")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>>())
      .def("run", &MultipoleMomentTask::run)
      .def_readwrite("settings", &MultipoleMomentTask::settings)
      .def_readwrite("generalSettings", &MultipoleMomentTask::generalSettings);
}
