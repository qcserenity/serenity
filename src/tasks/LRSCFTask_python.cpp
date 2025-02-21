/**
 * @file LRSCFTask_python.cpp
 *
 * @date Nov 03, 2020
 * @author Niklas Niemeyer
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
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

void export_LRSCFTask(py::module& spy) {
  py::class_<LRSCFTaskSettings>(spy, "LRSCFTaskSettings", "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("nEigen", &LRSCFTaskSettings::nEigen)
      .def_readwrite("conv", &LRSCFTaskSettings::conv)
      .def_readwrite("maxCycles", &LRSCFTaskSettings::maxCycles)
      .def_readwrite("maxSubspaceDimension", &LRSCFTaskSettings::maxSubspaceDimension)
      .def_readwrite("dominantThresh", &LRSCFTaskSettings::dominantThresh)
      .def_readwrite("func", &LRSCFTaskSettings::func)
      .def_readwrite("analysis", &LRSCFTaskSettings::analysis)
      .def_readwrite("besleyAtoms", &LRSCFTaskSettings::besleyAtoms)
      .def_readwrite("besleyCutoff", &LRSCFTaskSettings::besleyCutoff)
      .def_readwrite("excludeProjection", &LRSCFTaskSettings::excludeProjection)
      .def_readwrite("uncoupledSubspace", &LRSCFTaskSettings::uncoupledSubspace)
      .def_readwrite("fullFDEc", &LRSCFTaskSettings::fullFDEc)
      .def_readwrite("loadType", &LRSCFTaskSettings::loadType)
      .def_readwrite("gauge", &LRSCFTaskSettings::gauge)
      .def_readwrite("gaugeOrigin", &LRSCFTaskSettings::gaugeOrigin)
      .def_readwrite("frequencies", &LRSCFTaskSettings::frequencies)
      .def_readwrite("frequencyRange", &LRSCFTaskSettings::frequencyRange)
      .def_readwrite("damping", &LRSCFTaskSettings::damping)
      .def_readwrite("couplingPattern", &LRSCFTaskSettings::couplingPattern)
      .def_readwrite("method", &LRSCFTaskSettings::method)
      .def_readwrite("diis", &LRSCFTaskSettings::diis)
      .def_readwrite("diisStore", &LRSCFTaskSettings::diisStore)
      .def_readwrite("preopt", &LRSCFTaskSettings::preopt)
      .def_readwrite("cctrdens", &LRSCFTaskSettings::cctrdens)
      .def_readwrite("ccexdens", &LRSCFTaskSettings::ccexdens)
      .def_readwrite("sss", &LRSCFTaskSettings::sss)
      .def_readwrite("oss", &LRSCFTaskSettings::oss)
      .def_readwrite("nafThresh", &LRSCFTaskSettings::nafThresh)
      .def_readwrite("samedensity", &LRSCFTaskSettings::samedensity)
      .def_readwrite("subsystemgrid", &LRSCFTaskSettings::subsystemgrid)
      .def_readwrite("rpaScreening", &LRSCFTaskSettings::rpaScreening)
      .def_readwrite("restart", &LRSCFTaskSettings::restart)
      .def_readwrite("densFitJ", &LRSCFTaskSettings::densFitJ)
      .def_readwrite("densFitK", &LRSCFTaskSettings::densFitK)
      .def_readwrite("densFitLRK", &LRSCFTaskSettings::densFitLRK)
      .def_readwrite("densFitCache", &LRSCFTaskSettings::densFitCache)
      .def_readwrite("transitionCharges", &LRSCFTaskSettings::transitionCharges)
      .def_readwrite("partialResponseConstruction", &LRSCFTaskSettings::partialResponseConstruction)
      .def_readwrite("grimme", &LRSCFTaskSettings::grimme)
      .def_readwrite("adaptivePrescreening", &LRSCFTaskSettings::adaptivePrescreening)
      .def_readwrite("frozenCore", &LRSCFTaskSettings::frozenCore)
      .def_readwrite("frozenVirtual", &LRSCFTaskSettings::frozenVirtual)
      .def_readwrite("coreOnly", &LRSCFTaskSettings::coreOnly)
      .def_readwrite("ltconv", &LRSCFTaskSettings::ltconv)
      .def_readwrite("aocache", &LRSCFTaskSettings::aocache)
      .def_readwrite("triplet", &LRSCFTaskSettings::triplet)
      .def_readwrite("noCoupling", &LRSCFTaskSettings::noCoupling)
      .def_readwrite("approxCoulomb", &LRSCFTaskSettings::approxCoulomb)
      .def_readwrite("scfstab", &LRSCFTaskSettings::scfstab)
      .def_readwrite("stabroot", &LRSCFTaskSettings::stabroot)
      .def_readwrite("stabscal", &LRSCFTaskSettings::stabscal)
      .def_readwrite("noKernel", &LRSCFTaskSettings::noKernel)
      .def_readwrite("excGradList", &LRSCFTaskSettings::excGradList)
      .def_readwrite("hypthresh", &LRSCFTaskSettings::hypthresh)
      .def_readwrite("embedding", &LRSCFTaskSettings::embedding)
      .def_readwrite("grid", &LRSCFTaskSettings::grid)
      .def_readwrite("customFunc", &LRSCFTaskSettings::customFunc);

  py::class_<LRSCFTask<RESTRICTED>, std::shared_ptr<LRSCFTask<RESTRICTED>>>(spy, "LRSCFTask_R")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &LRSCFTask<RESTRICTED>::run)
      .def("getTransitions", &LRSCFTask<RESTRICTED>::getTransitions)
      .def("getProperties", &LRSCFTask<RESTRICTED>::getProperties)
      .def_readwrite("settings", &LRSCFTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &LRSCFTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<LRSCFTask<UNRESTRICTED>, std::shared_ptr<LRSCFTask<UNRESTRICTED>>>(spy, "LRSCFTask_U")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&, std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &LRSCFTask<UNRESTRICTED>::run)
      .def("getTransitions", &LRSCFTask<UNRESTRICTED>::getTransitions)
      .def("getProperties", &LRSCFTask<UNRESTRICTED>::getProperties)
      .def_readwrite("settings", &LRSCFTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &LRSCFTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
