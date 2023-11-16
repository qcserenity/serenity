/**
 * @file Libint_python.cpp
 *
 * @date Feb 26, 2018
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
#include "basis/BasisController.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

std::shared_ptr<Libint> libintConstr() {
  return Libint::getSharedPtr();
}

Eigen::MatrixXd comp2c1b1(std::shared_ptr<Libint> libint, LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis,
                          std::shared_ptr<Geometry> geo) {
  return libint->compute1eInts(op, basis, geo->getAtoms());
}

Eigen::MatrixXd comp2c1b2(std::shared_ptr<Libint> libint, LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis) {
  return libint->compute1eInts(op, basis);
}

Eigen::MatrixXd comp2c2b1(std::shared_ptr<Libint> libint, LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                          std::shared_ptr<BasisController> basis2, std::shared_ptr<Geometry> geo) {
  return libint->compute1eInts(op, basis1, basis2, geo->getAtoms());
}
Eigen::MatrixXd comp2c2b2(std::shared_ptr<Libint> libint, LIBINT_OPERATOR op, std::shared_ptr<BasisController> basis1,
                          std::shared_ptr<BasisController> basis2) {
  return libint->compute1eInts(op, basis1, basis2);
}

void export_Libint(py::module& spy) {
  py::enum_<LIBINT_OPERATOR>(spy, "INT_OPERATORS")
      .value("COULOMB", LIBINT_OPERATOR::coulomb)
      .value("NUCLEAR", LIBINT_OPERATOR::nuclear)
      .value("KINETIC", LIBINT_OPERATOR::kinetic)
      .value("OVERLAP", LIBINT_OPERATOR::overlap)
      .export_values();

  py::class_<Libint, std::shared_ptr<Libint>>(spy, "Libint", "Libint wrapped for basic integral evaluation")
      .def(py::init(&libintConstr))
      .def("compute", &comp2c1b1)
      .def("compute", &comp2c1b2)
      .def("compute", &comp2c2b1)
      .def("compute", &comp2c2b2);
}
