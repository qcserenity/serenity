/**
 * @file FXDTask_python.cpp
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
#include "tasks/FXDTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_FXDTask(py::module& spy) {
  py::class_<FXDTaskSettings>(spy, "FXDTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("loadType", &FXDTaskSettings::loadType)
      .def_readwrite("donoratoms", &FXDTaskSettings::donoratoms)
      .def_readwrite("acceptoratoms", &FXDTaskSettings::acceptoratoms)
      .def_readwrite("FED", &FXDTaskSettings::FED)
      .def_readwrite("FCD", &FXDTaskSettings::FCD)
      .def_readwrite("multistateFXD", &FXDTaskSettings::multistateFXD)
      .def_readwrite("states", &FXDTaskSettings::states)
      .def_readwrite("loewdinpopulation", &FXDTaskSettings::loewdinpopulation)
      .def_readwrite("writeTransformedExcitationVectors", &FXDTaskSettings::writeTransformedExcitationVectors);

  py::class_<FXDTask<Options::SCF_MODES::RESTRICTED>>(spy, "FXDTask_R")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &FXDTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &FXDTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &FXDTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<FXDTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "FXDTask_U")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}