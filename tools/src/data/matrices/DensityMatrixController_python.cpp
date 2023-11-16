/**
 * @file DensityMatrixController_python.cpp
 *
 * @date Mar 20, 2017
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
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
#include "energies/EnergyComponentController.h"
#include "grid/GridController.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd totalOnGrid(std::shared_ptr<DensityMatrixController<SCFMode>> dmatC,
                            std::shared_ptr<GridController> grid, double avgthrhld = 1e-11, double radialthrhld = 1e-9) {
  // generate data on grid
  auto bfOnGrid = std::make_shared<BasisFunctionOnGridController>(dmatC->getDensityMatrix().getBasisController(), grid,
                                                                  128, radialthrhld, 0);
  auto densOnGridCalc = std::make_shared<DensityOnGridCalculator<SCFMode>>(bfOnGrid, avgthrhld);
  DensityMatrixDensityOnGridController<SCFMode> dmdogCont(densOnGridCalc, dmatC, 0);

  auto& densOnGrid = dmdogCont.getDensityOnGrid();
  Eigen::VectorXd list(grid->getNGridPoints());
  for (unsigned int i = 0; i < grid->getNGridPoints(); ++i) {
    double tmp = 0.0;
    for_spin(densOnGrid) {
      tmp += densOnGrid_spin[i];
    };
    list[i] = tmp;
  }
  return list;
}

Eigen::MatrixXd totalU(std::shared_ptr<DensityMatrixController<UNRESTRICTED>> dmatC) {
  return dmatC->getDensityMatrix().total();
}
Eigen::MatrixXd alphaU(std::shared_ptr<DensityMatrixController<UNRESTRICTED>> dmatC) {
  return dmatC->getDensityMatrix().alpha;
}
Eigen::MatrixXd betaU(std::shared_ptr<DensityMatrixController<UNRESTRICTED>> dmatC) {
  return dmatC->getDensityMatrix().beta;
}

void setR(std::shared_ptr<DensityMatrixController<RESTRICTED>> dmatC, Eigen::MatrixXd& newMat) {
  auto mat = dmatC->getDensityMatrix();
  mat = newMat;
  dmatC->setDensityMatrix(mat);
}

void export_DensityMatrixController(py::module& spy) {
  py::class_<DensityMatrixController<Options::SCF_MODES::RESTRICTED>, std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>>(
      spy, "DensityMatrix_R", "@brief A restricted density matrix.")
      .def("totalOnGrid", totalOnGrid<Options::SCF_MODES::RESTRICTED>, py::arg("grid"), py::arg("avgthrhld") = 1e-11,
           py::arg("radialthrhld") = 1e-9)
      .def("total", &DensityMatrixController<Options::SCF_MODES::RESTRICTED>::getDensityMatrix)
      .def("set", &setR);
  py::class_<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>,
             std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>>(
      spy, "DensityMatrix_U", "@brief An unrestricted density matrix.")
      .def("totalOnGrid", totalOnGrid<Options::SCF_MODES::UNRESTRICTED>, py::arg("grid"), py::arg("avgthrhld") = 1e-11,
           py::arg("radialthrhld") = 1e-9)
      .def("total", &totalU)
      .def("alpha", &alphaU)
      .def("beta", &alphaU);
}