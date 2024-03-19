/**
 * @file GeometryOptimizationTask_python.cpp
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
#include "tasks/GeometryOptimizationTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_GeometryOptimizationTask(py::module& spy) {
  py::class_<GeometryOptimizationTaskSettings>(
      spy, "GeometryOptimizationTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("gradType", &GeometryOptimizationTaskSettings::gradType)
      .def_readwrite("maxCycles", &GeometryOptimizationTaskSettings::maxCycles)
      .def_readwrite("rmsgradThresh", &GeometryOptimizationTaskSettings::rmsgradThresh)
      .def_readwrite("energyChangeThresh", &GeometryOptimizationTaskSettings::energyChangeThresh)
      .def_readwrite("maxGradThresh", &GeometryOptimizationTaskSettings::maxGradThresh)
      .def_readwrite("stepThresh", &GeometryOptimizationTaskSettings::stepThresh)
      .def_readwrite("maxStepThresh", &GeometryOptimizationTaskSettings::maxStepThresh)
      .def_readwrite("numGradStepSize", &GeometryOptimizationTaskSettings::numGradStepSize)
      .def_readwrite("transInvar", &GeometryOptimizationTaskSettings::transInvar)
      .def_readwrite("FaTmaxCycles", &GeometryOptimizationTaskSettings::FaTmaxCycles)
      .def_readwrite("FaTenergyConvThresh", &GeometryOptimizationTaskSettings::FaTConvThresh)
      .def_readwrite("FaTgridCutOff", &GeometryOptimizationTaskSettings::FaTgridCutOff)
      .def_readwrite("embedding", &GeometryOptimizationTaskSettings::embedding);

  py::class_<GeometryOptimizationTask<RESTRICTED>>(spy, "GeometryOptimizationTask_R")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GeometryOptimizationTask<RESTRICTED>::run)
      .def_readwrite("settings", &GeometryOptimizationTask<RESTRICTED>::settings)
      .def_readwrite("generalSettings", &GeometryOptimizationTask<RESTRICTED>::generalSettings);
  py::class_<GeometryOptimizationTask<UNRESTRICTED>>(spy, "GeometryOptimizationTask_U")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GeometryOptimizationTask<UNRESTRICTED>::run)
      .def_readwrite("settings", &GeometryOptimizationTask<UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &GeometryOptimizationTask<UNRESTRICTED>::generalSettings);
}