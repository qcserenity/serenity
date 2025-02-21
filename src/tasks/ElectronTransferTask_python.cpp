/**
 * @file ElectronTransferTask_python.cpp
 *
 * @date Nov 11, 2024
 * @author: Lukas Lampe
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
#include "postHF/ET/FDEETCalculator.h"
#include "system/SystemController.h"
#include "tasks/ElectronTransferTask.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

double coupling(ElectronTransferTask etTask) {
  return etTask.getFDEETCalculator()->getAnalyticalCoupling();
}

Eigen::VectorXd eigenvalues(ElectronTransferTask etTask) {
  return etTask.getFDEETCalculator()->getEigenValues();
}

Eigen::MatrixXd hamiltonian(ElectronTransferTask etTask) {
  return etTask.getFDEETCalculator()->getHamiltonian();
}

void export_ElectronTransferTask(py::module& spy) {
  py::class_<ElectronTransferTaskSettings>(spy, "ElectronTransferTaskSettings",
                                           "@brief Default constructor for Settings all set to their default values.")
      .def_readwrite("couple", &ElectronTransferTaskSettings::couple)
      .def_readwrite("states", &ElectronTransferTaskSettings::states)
      .def_readwrite("population", &ElectronTransferTaskSettings::population)
      .def_readwrite("spindensity", &ElectronTransferTaskSettings::spindensity)
      .def_readwrite("spinpopulation", &ElectronTransferTaskSettings::spinpopulation)
      .def_readwrite("disjoint", &ElectronTransferTaskSettings::disjoint)
      .def_readwrite("diskMode", &ElectronTransferTaskSettings::diskMode)
      .def_readwrite("printContributions", &ElectronTransferTaskSettings::printContributions)
      .def_readwrite("useHFCoupling", &ElectronTransferTaskSettings::useHFCoupling)
      .def_readwrite("coupleAdiabaticStates", &ElectronTransferTaskSettings::coupleAdiabaticStates);

  py::class_<ElectronTransferTask, std::shared_ptr<ElectronTransferTask>>(spy, "ElectronTransferTask")
      .def(py::init<std::vector<std::shared_ptr<SystemController>>&>())
      .def("run", &ElectronTransferTask::run)
      .def("coupling", coupling)
      .def("eigenvalues", eigenvalues)
      .def("hamiltonian", hamiltonian)
      .def_readwrite("settings", &ElectronTransferTask::settings)
      .def_readwrite("generalSettings", &ElectronTransferTask::generalSettings);
}