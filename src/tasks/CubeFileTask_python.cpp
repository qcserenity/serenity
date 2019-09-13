/**
 * @file   CubeFileTask_python.cpp
 *
 * @date   Sep 21, 2016
 * @author Jan Unsleber
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
#include "tasks/CubeFileTask.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_CubeFileTask(py::module &spy){

  py::class_<CubeFileTaskSettings>(spy,"CubeFileTaskSettings",
    "@brief Default constructor for Settings all set to their default values.")
  .def_readwrite("density",&CubeFileTaskSettings::density)
  .def_readwrite("allOrbitals",&CubeFileTaskSettings::allOrbitals)
  .def_readwrite("cubeSpacing",&CubeFileTaskSettings::cubeSpacing)
  .def_readwrite("cubeBorder",&CubeFileTaskSettings::cubeBorder)
  .def_readwrite("orbitals",&CubeFileTaskSettings::orbitals)
  .def_readwrite("electrostaticPot",&CubeFileTaskSettings::electrostaticPot)
  .def_readwrite("sedd",&CubeFileTaskSettings::sedd)
  .def_readwrite("dori",&CubeFileTaskSettings::dori)
  .def_readwrite("signedDensity",&CubeFileTaskSettings::signedDensity)
  .def_readwrite("allOrbitals",&CubeFileTaskSettings::allOrbitals)
  .def_readwrite("occOrbitals",&CubeFileTaskSettings::occOrbitals)
  .def_readwrite("electrostaticPot",&CubeFileTaskSettings::electrostaticPot);

  py::class_<CubeFileTask<RESTRICTED> >(spy,"CubeFileTask_R")
        .def(py::init<const std::vector<std::shared_ptr<SystemController> > &,
                      const std::vector<std::shared_ptr<SystemController> > &>())
        .def("run", &CubeFileTask<RESTRICTED>::run )
        .def_readwrite("settings",&CubeFileTask<Options::SCF_MODES::RESTRICTED>::settings);
  py::class_<CubeFileTask<UNRESTRICTED> >(spy,"CubeFileTask_U")
        .def(py::init<const std::vector<std::shared_ptr<SystemController> > &,
                      const std::vector<std::shared_ptr<SystemController> > &>())
        .def("run", &CubeFileTask<UNRESTRICTED>::run )
        .def_readwrite("settings",&CubeFileTask<Options::SCF_MODES::UNRESTRICTED>::settings);
}






