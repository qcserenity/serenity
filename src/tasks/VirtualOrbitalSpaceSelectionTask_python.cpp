/**
 * @file VirtualOrbitalSpaceSelectionTask_python.cpp
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
#include "tasks/VirtualOrbitalSpaceSelectionTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_VirtualOrbitalSpaceSelectionTask(py::module& spy) {
  py::class_<VirtualOrbitalSpaceSelectionTaskSettings>(
      spy, "VirtualOrbitalSpaceSelectionTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("excludeProjection", &VirtualOrbitalSpaceSelectionTaskSettings::excludeProjection)
      .def_readwrite("localCanonicalVirtuals", &VirtualOrbitalSpaceSelectionTaskSettings::localCanonicalVirtuals)
      .def_readwrite("envCanonicalVirtuals", &VirtualOrbitalSpaceSelectionTaskSettings::envCanonicalVirtuals)
      .def_readwrite("localizedVirtualorbitals", &VirtualOrbitalSpaceSelectionTaskSettings::localizedVirtualorbitals)
      .def_readwrite("localizedEnvVirtualorbitals", &VirtualOrbitalSpaceSelectionTaskSettings::localizedEnvVirtualorbitals)
      .def_readwrite("recalculateFockMatrix", &VirtualOrbitalSpaceSelectionTaskSettings::recalculateFockMatrix)
      .def_readwrite("identifier", &VirtualOrbitalSpaceSelectionTaskSettings::identifier)
      .def_readwrite("mixingOccAndVirtualOrbitals", &VirtualOrbitalSpaceSelectionTaskSettings::mixingOccAndVirtualOrbitals)
      .def_readwrite("relaxation", &VirtualOrbitalSpaceSelectionTaskSettings::relaxation)
      .def_readwrite("onlyOne", &VirtualOrbitalSpaceSelectionTaskSettings::onlyOne)
      .def_readwrite("embedding", &VirtualOrbitalSpaceSelectionTaskSettings::embedding);

  py::class_<VirtualOrbitalSpaceSelectionTask<RESTRICTED>, std::shared_ptr<VirtualOrbitalSpaceSelectionTask<RESTRICTED>>>(
      spy, "VirtualOrbitalSpaceSelectionTask_R")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &VirtualOrbitalSpaceSelectionTask<RESTRICTED>::run)
      .def_readwrite("settings", &VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<VirtualOrbitalSpaceSelectionTask<UNRESTRICTED>, std::shared_ptr<VirtualOrbitalSpaceSelectionTask<UNRESTRICTED>>>(
      spy, "VirtualOrbitalSpaceSelectionTask_U")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &VirtualOrbitalSpaceSelectionTask<UNRESTRICTED>::run)
      .def_readwrite("settings", &VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &VirtualOrbitalSpaceSelectionTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}