/**
 * @file GWTask_python.cpp
 *
 * @date Feb 14, 2025
 * @author Johannes Tölle
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
#include "tasks/GWTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_GWTask(py::module& spy) {
  py::class_<GWTaskSettings>(spy, "GWTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("mbpttype", &GWTaskSettings::mbpttype)
      .def_readwrite("gwtype", &GWTaskSettings::gwtype)
      .def_readwrite("linearized", &GWTaskSettings::linearized)
      .def_readwrite("qpiterations", &GWTaskSettings::qpiterations)
      .def_readwrite("eta", &GWTaskSettings::eta)
      .def_readwrite("nVirt", &GWTaskSettings::nVirt)
      .def_readwrite("nOcc", &GWTaskSettings::nOcc)
      .def_readwrite("integrationPoints", &GWTaskSettings::integrationPoints)
      .def_readwrite("padePoints", &GWTaskSettings::padePoints)
      .def_readwrite("fermiShift", &GWTaskSettings::fermiShift)
      .def_readwrite("derivativeShift", &GWTaskSettings::derivativeShift)
      .def_readwrite("imagShift", &GWTaskSettings::imagShift)
      .def_readwrite("gridCutOff", &GWTaskSettings::gridCutOff)
      .def_readwrite("evGW", &GWTaskSettings::evGW)
      .def_readwrite("evGWcycles", &GWTaskSettings::evGWcycles)
      .def_readwrite("diis", &GWTaskSettings::diis)
      .def_readwrite("diisMaxStore", &GWTaskSettings::diisMaxStore)
      .def_readwrite("ConvergenceThreshold", &GWTaskSettings::ConvergenceThreshold)
      .def_readwrite("nafThresh", &GWTaskSettings::nafThresh)
      .def_readwrite("subsystemAuxillaryBasisOnly", &GWTaskSettings::subsystemAuxillaryBasisOnly)
      .def_readwrite("hybrid", &GWTaskSettings::hybrid)
      .def_readwrite("freq", &GWTaskSettings::freq)
      .def_readwrite("damping", &GWTaskSettings::damping)
      .def_readwrite("gap", &GWTaskSettings::gap)
      .def_readwrite("environmentScreening", &GWTaskSettings::environmentScreening)
      .def_readwrite("ltconv", &GWTaskSettings::ltconv)
      .def_readwrite("frozenCore", &GWTaskSettings::frozenCore)
      .def_readwrite("coreOnly", &GWTaskSettings::coreOnly)
      .def_readwrite("densFitCache", &GWTaskSettings::densFitCache);

  py::class_<GWTask<RESTRICTED>, std::shared_ptr<GWTask<RESTRICTED>>>(spy, "GWTask_R")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GWTask<RESTRICTED>::run)
      .def_readwrite("settings", &GWTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &GWTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<GWTask<UNRESTRICTED>, std::shared_ptr<GWTask<UNRESTRICTED>>>(spy, "GWTask_U")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &GWTask<UNRESTRICTED>::run)
      .def_readwrite("settings", &GWTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &GWTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
