/**
 * @file GridController_python.cpp
 *
 * @date Mar 19, 2017
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
#include "geometry/Point.h"
#include "grid/Grid.h"
#include "grid/GridController.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Serenity;

static std::shared_ptr<GridController> createGridControllerPtr1(std::vector<double> x, std::vector<double> y,
                                                                std::vector<double> z, std::vector<double> w) {
  assert(x.size() == y.size());
  assert(x.size() == z.size());
  assert(x.size() == w.size());
  std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, x.size()));
  for (unsigned int i = 0; i < (unsigned int)x.size(); ++i) {
    (*points)(0, i) = x[i];
    (*points)(1, i) = y[i];
    (*points)(2, i) = z[i];
  }
  Eigen::VectorXd tmp2 = Eigen::Map<Eigen::VectorXd>(w.data(), w.size());
  std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(tmp2));
  return std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
}
static std::shared_ptr<GridController> createGridControllerPtr2(std::vector<double> x, std::vector<double> y,
                                                                std::vector<double> z) {
  assert(x.size() == y.size());
  assert(x.size() == z.size());
  std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, x.size()));
  for (unsigned int i = 0; i < (unsigned int)x.size(); ++i) {
    (*points)(0, i) = x[i];
    (*points)(1, i) = y[i];
    (*points)(2, i) = z[i];
  }
  std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(x.size()));
  weights->setOnes();
  return std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
}

void export_GridController(py::module& spy) {
  py::class_<GridController, std::shared_ptr<GridController>>(spy, "Grid")
      .def(py::init(&createGridControllerPtr1))
      .def(py::init(&createGridControllerPtr2))
      .def("getGridPoints", &GridController::getGridPoints)
      .def("getWeights", &GridController::getWeights);
}