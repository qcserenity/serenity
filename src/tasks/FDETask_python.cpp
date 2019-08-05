/**
 * @file FDETask_python.cpp
 *
 * @date Apr 28, 2016
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
#include "tasks/FDETask.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_FDETask(py::module &spy){

  py::class_<FDETaskSettings>(spy,"FDETaskSettings",
      "@brief Default constructor for Settings all set to their default values.")
    .def_readwrite("naddKinFunc",&FDETaskSettings::naddKinFunc)
    .def_readwrite("naddXCFunc",&FDETaskSettings::naddXCFunc)
    .def_readwrite("embeddingMode",&FDETaskSettings::embeddingMode)
    .def_readwrite("gridCutOff",&FDETaskSettings::gridCutOff)
    .def_readwrite("naddDispersion",&FDETaskSettings::dispersion);

  py::class_<FDETask<Options::SCF_MODES::RESTRICTED> >(spy,"FDETask_R")
      .def(py::init<std::shared_ptr<SystemController>,
                    std::vector<std::shared_ptr<SystemController> > &>())
      .def("run", &FDETask<Options::SCF_MODES::RESTRICTED>::run )
      .def_readwrite("settings",&FDETask<Options::SCF_MODES::RESTRICTED>::settings);
  py::class_<FDETask<Options::SCF_MODES::UNRESTRICTED> >(spy,"FDETask_U")
      .def(py::init<std::shared_ptr<SystemController>,
                    std::vector<std::shared_ptr<SystemController> > &>())
      .def("run", &FDETask<Options::SCF_MODES::UNRESTRICTED>::run )
      .def_readwrite("settings",&FDETask<Options::SCF_MODES::UNRESTRICTED>::settings);
}


