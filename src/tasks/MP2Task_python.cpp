/**
 * @file MP2Task_python.cpp
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
#include "tasks/MP2Task.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_MP2Task(py::module& spy) {
  py::class_<MP2TaskSettings>(spy, "MP2TaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("mp2Type", &MP2TaskSettings::mp2Type)
      .def_readwrite("sss", &MP2TaskSettings::sss)
      .def_readwrite("oss", &MP2TaskSettings::oss)
      .def_readwrite("ltconv", &MP2TaskSettings::ltconv)
      .def_readwrite("maxResidual", &MP2TaskSettings::maxResidual)
      .def_readwrite("maxCycles", &MP2TaskSettings::maxCycles)
      .def_readwrite("writePairEnergies", &MP2TaskSettings::writePairEnergies)
      .def_readwrite("unrelaxedDensity", &MP2TaskSettings::unrelaxedDensity)
      .def_readwrite("lcSettings", &MP2TaskSettings::lcSettings);

  py::class_<MP2Task<Options::SCF_MODES::RESTRICTED>>(spy, "MP2Task_R")
      .def(py::init<std::shared_ptr<SystemController>, std::vector<std::shared_ptr<SystemController>>&>())
      .def_readwrite("settings", &MP2Task<Options::SCF_MODES::RESTRICTED>::settings)
      .def("run", &MP2Task<Options::SCF_MODES::RESTRICTED>::run)
      .def_readwrite("generalSettings", &MP2Task<Options::SCF_MODES::RESTRICTED>::generalSettings);
}