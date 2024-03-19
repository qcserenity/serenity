/**
 * @file GradientTask_python.cpp
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
#include "tasks/GradientTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_GradientTask(py::module& spy) {
  py::class_<GradientTaskSettings>(spy, "GradientTaskSettings",
                                   "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("gradType", &GradientTaskSettings::gradType)
      .def_readwrite("numGradStepSize", &GradientTaskSettings::numGradStepSize)
      .def_readwrite("transInvar", &GradientTaskSettings::transInvar)
      .def_readwrite("fatMaxCycles", &GradientTaskSettings::FaTmaxCycles)
      .def_readwrite("fatEnergyConvThresh", &GradientTaskSettings::FaTenergyConvThresh)
      .def_readwrite("fdeGridCutOff", &GradientTaskSettings::FDEgridCutOff)
      .def_readwrite("print", &GradientTaskSettings::print)
      .def_readwrite("printTotal", &GradientTaskSettings::printTotal)
      .def_readwrite("embedding", &GradientTaskSettings::embedding);

  py::class_<GradientTask<Options::SCF_MODES::RESTRICTED>>(spy, "GradientTask_R")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GradientTask<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("settings", &GradientTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &GradientTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<GradientTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "GradientTask_U")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GradientTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def_readwrite("settings", &GradientTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &GradientTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
