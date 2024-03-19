/**
 * @file CoupledClusterTask_python.cpp
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
#include "tasks/CoupledClusterTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;
void export_CoupledClusterTask(py::module& spy) {
  py::class_<CoupledClusterTaskSettings>(spy, "CoupledClusterTaskSettings",
                                         "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("level", &CoupledClusterTaskSettings::level)
      .def_readwrite("maxCycles", &CoupledClusterTaskSettings::maxCycles)
      .def_readwrite("normThreshold", &CoupledClusterTaskSettings::normThreshold)
      .def_readwrite("writePairEnergies", &CoupledClusterTaskSettings::writePairEnergies);

  py::class_<CoupledClusterTask>(spy, "CoupledClusterTask_R")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &CoupledClusterTask::run);
  //    class_<CoupledClusterTask<Options::SCF_MODES::UNRESTRICTED> >("CoupledClusterTask_U",
  //    init<std::shared_ptr<SystemController> >())
  //        .def("run", &CoupledClusterTask<Options::SCF_MODES::RESTRICTED>::run )
  //        ;
}