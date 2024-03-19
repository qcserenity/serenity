/**
 * @file GeneralizedDOSTask_python.cpp
 *
 * @date Sep 07, 2022
 * @author: Moritz Bensberg
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
#include "analysis/directOrbitalSelection/DirectOrbitalSelection.h"
#include "system/SystemController.h"
#include "tasks/GeneralizedDOSTask.h"
/* Include Std and External Headers */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

std::vector<std::vector<std::vector<unsigned int>>> getOrbitalGroupIndices(GeneralizedDOSTask<RESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getOrbitalGroups();
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

std::vector<std::vector<std::vector<unsigned int>>>
getOrbitalGroupIndicesAlpha(GeneralizedDOSTask<UNRESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getOrbitalGroups().alpha;
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

std::vector<std::vector<std::vector<unsigned int>>> getOrbitalGroupIndicesBeta(GeneralizedDOSTask<UNRESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getOrbitalGroups().beta;
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

std::vector<std::vector<std::vector<unsigned int>>>
getUnmappableOrbitalGroupIndices(GeneralizedDOSTask<RESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getUnmappableOrbitalGroups();
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

std::vector<std::vector<std::vector<unsigned int>>>
getUnmappableOrbitalGroupIndicesAlpha(GeneralizedDOSTask<UNRESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getUnmappableOrbitalGroups().alpha;
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

std::vector<std::vector<std::vector<unsigned int>>>
getUnmappableOrbitalGroupIndicesBeta(GeneralizedDOSTask<UNRESTRICTED>& generalizedDosTask) {
  const auto& groups = generalizedDosTask.getUnmappableOrbitalGroups().beta;
  return DOSOrbitalGroup::getIndicesFromGroups(groups);
}

void export_GeneralizedDOSTask(py::module& spy) {
  py::class_<GeneralizedDOSTaskSettings>(spy, "GeneralizedDOSTaskSettings",
                                         "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("similarityLocThreshold", &GeneralizedDOSTaskSettings::similarityLocThreshold)
      .def_readwrite("similarityKinEnergyThreshold", &GeneralizedDOSTaskSettings::similarityKinEnergyThreshold)
      .def_readwrite("prioFirst", &GeneralizedDOSTaskSettings::prioFirst)
      .def_readwrite("localizationThreshold", &GeneralizedDOSTaskSettings::localizationThreshold)
      .def_readwrite("populationAlgorithm", &GeneralizedDOSTaskSettings::populationAlgorithm)
      .def_readwrite("checkDegeneracies", &GeneralizedDOSTaskSettings::checkDegeneracies)
      .def_readwrite("degeneracyFactor", &GeneralizedDOSTaskSettings::degeneracyFactor)
      .def_readwrite("usePiBias", &GeneralizedDOSTaskSettings::usePiBias)
      .def_readwrite("biasThreshold", &GeneralizedDOSTaskSettings::biasThreshold)
      .def_readwrite("biasAverage", &GeneralizedDOSTaskSettings::biasAverage)
      .def_readwrite("writeScores", &GeneralizedDOSTaskSettings::writeScores)
      .def_readwrite("scoreStart", &GeneralizedDOSTaskSettings::scoreStart)
      .def_readwrite("scoreEnd", &GeneralizedDOSTaskSettings::scoreEnd)
      .def_readwrite("nTest", &GeneralizedDOSTaskSettings::nTest)
      .def_readwrite("mapVirtuals", &GeneralizedDOSTaskSettings::mapVirtuals)
      .def_readwrite("writeGroupsToFile", &GeneralizedDOSTaskSettings::writeGroupsToFile)
      .def_readwrite("bestMatchMapping", &GeneralizedDOSTaskSettings::bestMatchMapping);

  py::class_<GeneralizedDOSTask<Options::SCF_MODES::RESTRICTED>>(spy, "GeneralizedDOSTask_R")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>, std::vector<std::shared_ptr<SystemController>>>())
      .def("run", &GeneralizedDOSTask<Options::SCF_MODES::RESTRICTED>::run)
      .def("getOrbitalGroupIndices", &getOrbitalGroupIndices)
      .def("getUnmappableOrbitalGroupIndices", &getOrbitalGroupIndices)
      .def_readwrite("settings", &GeneralizedDOSTask<Options::SCF_MODES::RESTRICTED>::settings)
      .def_readwrite("generalSettings", &GeneralizedDOSTask<Options::SCF_MODES::RESTRICTED>::generalSettings);
  py::class_<GeneralizedDOSTask<Options::SCF_MODES::UNRESTRICTED>>(spy, "GeneralizedDOSTask_U")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>, std::vector<std::shared_ptr<SystemController>>>())
      .def("run", &GeneralizedDOSTask<Options::SCF_MODES::UNRESTRICTED>::run)
      .def("getOrbitalGroupIndicesAlpha", &getOrbitalGroupIndicesAlpha)
      .def("getOrbitalGroupIndicesBeta", &getOrbitalGroupIndicesBeta)
      .def("getUnmappableOrbitalGroupIndicesAlpha", &getOrbitalGroupIndicesAlpha)
      .def("getUnmappableOrbitalGroupIndicesBeta", &getOrbitalGroupIndicesBeta)
      .def_readwrite("settings", &GeneralizedDOSTask<Options::SCF_MODES::UNRESTRICTED>::settings)
      .def_readwrite("generalSettings", &GeneralizedDOSTask<Options::SCF_MODES::UNRESTRICTED>::generalSettings);
}
