/**
 * @file   PlotTask_python.cpp
 *
 * @date   Sep 21, 2016
 * @author Jan Unsleber
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
#include "tasks/PlotTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_PlotTask(py::module& spy) {
  py::class_<PlotTaskSettings>(spy, "PlotTaskSettings",
                               "@brief Default constructor for Settings all "
                               "set to their default values.")

      .def_readwrite("p1", &PlotTaskSettings::p1)
      .def_readwrite("p2", &PlotTaskSettings::p2)
      .def_readwrite("p3", &PlotTaskSettings::p3)
      .def_readwrite("p4", &PlotTaskSettings::p4)
      .def_readwrite("atom1", &PlotTaskSettings::atom1)
      .def_readwrite("atom2", &PlotTaskSettings::atom2)
      .def_readwrite("atom3", &PlotTaskSettings::atom3)
      .def_readwrite("atom4", &PlotTaskSettings::atom4)
      .def_readwrite("gridSpacing", &PlotTaskSettings::gridSpacing)
      .def_readwrite("borderWidth", &PlotTaskSettings::borderWidth)
      .def_readwrite("xUnitVector", &PlotTaskSettings::xUnitVector)
      .def_readwrite("yUnitVector", &PlotTaskSettings::yUnitVector)
      .def_readwrite("zUnitVector", &PlotTaskSettings::zUnitVector)
      .def_readwrite("projectCutOffRadius", &PlotTaskSettings::projectCutOffRadius)
      .def_readwrite("xyHeatmap", &PlotTaskSettings::xyHeatmap)
      .def_readwrite("density", &PlotTaskSettings::density)
      .def_readwrite("allOrbitals", &PlotTaskSettings::allOrbitals)
      .def_readwrite("occOrbitals", &PlotTaskSettings::occOrbitals)
      .def_readwrite("electrostaticPot", &PlotTaskSettings::electrostaticPot)
      .def_readwrite("sedd", &PlotTaskSettings::sedd)
      .def_readwrite("dori", &PlotTaskSettings::dori)
      .def_readwrite("elf", &PlotTaskSettings::elf)
      .def_readwrite("elfts", &PlotTaskSettings::elfts)
      .def_readwrite("signedDensity", &PlotTaskSettings::signedDensity)
      .def_readwrite("orbitals", &PlotTaskSettings::orbitals)
      .def_readwrite("maxGridPoints", &PlotTaskSettings::maxGridPoints)
      .def_readwrite("cavity", &PlotTaskSettings::cavity)
      .def_readwrite("gridCoordinates", &PlotTaskSettings::gridCoordinates)
      .def_readwrite("ntos", &PlotTaskSettings::ntos)
      .def_readwrite("ntoPlotThreshold", &PlotTaskSettings::ntoPlotThreshold)
      .def_readwrite("excitations", &PlotTaskSettings::excitations)
      .def_readwrite("nros", &PlotTaskSettings::nros)
      .def_readwrite("nrominimum", &PlotTaskSettings::nrominimum)
      .def_readwrite("cctrdens", &PlotTaskSettings::cctrdens)
      .def_readwrite("ccexdens", &PlotTaskSettings::ccexdens);

  py::class_<PlotTask<RESTRICTED>>(spy, "PlotTask_R")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &PlotTask<RESTRICTED>::run)
      .def_readwrite("settings", &PlotTask<Options::SCF_MODES::RESTRICTED>::settings);
  py::class_<PlotTask<UNRESTRICTED>>(spy, "PlotTask_U")
      .def(py::init<const std::vector<std::shared_ptr<SystemController>>&, const std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &PlotTask<UNRESTRICTED>::run)
      .def_readwrite("settings", &PlotTask<Options::SCF_MODES::UNRESTRICTED>::settings);
}
