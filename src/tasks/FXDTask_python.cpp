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
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/algorithm/string.hpp>

namespace py = pybind11;
using namespace Serenity;

template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd& getMultistateExcitations(const FXDTask<SCFMode>& fxd, std::string type) {
  boost::algorithm::to_lower(type);
  if (type == "ctad")
    return fxd._ctAD;
  if (type == "ctda")
    return fxd._ctDA;
  if (type == "lea")
    return fxd._leA;
  if (type == "led")
    return fxd._leD;
  throw SerenityError("FXDTask: Invalid type for excitation!" + type);
}

template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd& getMultistateCouplings(const FXDTask<SCFMode>& fxd, std::string type) {
  boost::algorithm::to_lower(type);
  if (type == "ctadxctda")
    return fxd._ctADxctDA;
  if (type == "leaxled")
    return fxd._leAxleD;
  if (type == "ctadxlea")
    return fxd._ctADxleA;
  if (type == "ctdaxlea")
    return fxd._ctDAxleA;
  if (type == "ctadxled")
    return fxd._ctADxleD;
  if (type == "ctdaxled")
    return fxd._ctDAxleD;
  throw SerenityError("FXDTask: Invalid type for coupling!" + type);
}

template<Options::SCF_MODES SCFMode>
const std::map<std::pair<unsigned int, unsigned int>, Eigen::Vector4d>& getFEDResults(const FXDTask<SCFMode>& fxd) {
  return fxd._fedResults;
}

template<Options::SCF_MODES SCFMode>
const std::map<std::pair<unsigned int, unsigned int>, Eigen::Vector4d>& getFCDResults(const FXDTask<SCFMode>& fxd) {
  return fxd._fcdResults;
}

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
      .def("getMultistateExcitations", &getMultistateExcitations<RESTRICTED>)
      .def("getMultistateCouplings", &getMultistateCouplings<RESTRICTED>)
      .def("getFEDResults", &getFEDResults<RESTRICTED>,
           "@brief Get the FED results (after running an FXDTask_R with fed set to true)\n"
           "@returns A map with the 1-based excitation indices as key and the FED results as value\n"
           "The FED results are a vector of length 4 with the overall coupling as the fourth element.")
      .def("getFCDResults", &getFCDResults<RESTRICTED>)
      .def_readwrite("settings", &FXDTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &FXDTask<Options::SCF_MODES::RESTRICTED>::generalSettings);

  py::class_<FXDTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "FXDTask_U")
      .def(py::init<std::shared_ptr<SystemController>>())
      .def("run", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def("getMultistateExcitations", &getMultistateExcitations<UNRESTRICTED>)
      .def("getMultistateCouplings", &getMultistateCouplings<UNRESTRICTED>)
      .def("getFEDResults", &getFEDResults<UNRESTRICTED>,
           "@brief Get the FED results (after running an FXDTask_U with fed set to true)\n"
           "@returns A map with the 1-based excitation indices as key and the FED results as value\n"
           "The FED results are a vector of length 4 with the overall coupling as the fourth element.")
      .def("getFCDResults", &getFCDResults<UNRESTRICTED>)
      .def_readwrite("settings", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &FXDTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}