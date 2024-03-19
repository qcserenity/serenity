/**
 * @file FreezeAndThawTask_python.cpp
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
#include "tasks/FreezeAndThawTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_FreezeAndThawTask(py::module& spy) {
  py::class_<FreezeAndThawTaskSettings>(spy, "FreezeAndThawTaskSettings",
                                        "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("maxCycles", &FreezeAndThawTaskSettings::maxCycles)
      .def_readwrite("convThresh", &FreezeAndThawTaskSettings::convThresh)
      .def_readwrite("gridCutOff", &FreezeAndThawTaskSettings::gridCutOff)
      .def_readwrite("smallSupersystemGrid", &FreezeAndThawTaskSettings::smallSupersystemGrid)
      .def_readwrite("basisExtThresh", &FreezeAndThawTaskSettings::basisExtThresh)
      .def_readwrite("extendBasis", &FreezeAndThawTaskSettings::extendBasis)
      .def_readwrite("useConvAcceleration", &FreezeAndThawTaskSettings::useConvAcceleration)
      .def_readwrite("diisStart", &FreezeAndThawTaskSettings::diisStart)
      .def_readwrite("diisEnd", &FreezeAndThawTaskSettings::diisEnd)
      .def_readwrite("calculateSolvationEnergy", &FreezeAndThawTaskSettings::calculateSolvationEnergy)
      .def_readwrite("calculateUnrelaxedMP2Density", &FreezeAndThawTaskSettings::calculateUnrelaxedMP2Density)
      .def_readwrite("mp2Type", &FreezeAndThawTaskSettings::mp2Type)
      .def_readwrite("keepCoulombCache", &FreezeAndThawTaskSettings::keepCoulombCache)
      .def_readwrite("finalEnergyEvaluation", &FreezeAndThawTaskSettings::finalEnergyEvaluation)
      .def_readwrite("embedding", &FreezeAndThawTaskSettings::embedding)
      .def_readwrite("lcSettings", &FreezeAndThawTaskSettings::lcSettings);

  py::class_<FreezeAndThawTask<RESTRICTED>, std::shared_ptr<FreezeAndThawTask<RESTRICTED>>>(spy, "FreezeAndThawTask_R")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &FreezeAndThawTask<RESTRICTED>::run)
      .def_readwrite("settings", &FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &FreezeAndThawTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<FreezeAndThawTask<UNRESTRICTED>, std::shared_ptr<FreezeAndThawTask<UNRESTRICTED>>>(spy,
                                                                                                "FreezeAndThawTask_U")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &FreezeAndThawTask<UNRESTRICTED>::run)
      .def_readwrite("settings", &FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}