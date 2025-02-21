/**
 * @file Geometry_python.cpp
 *
 * @date Apr 25, 2016
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
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "geometry/XyzFileToGeometryConverter.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace Serenity;

std::shared_ptr<Geometry> createGeometryPtr(std::string path) {
  XyzFileToGeometryConverter converter(path);
  return std::shared_ptr<Geometry>(converter.readGeometry());
}

void export_Geometry(py::module& spy) {
  py::class_<Geometry, std::shared_ptr<Geometry>>(
      spy, "Geometry",
      "A geometry consisting of a list of atoms.\n"
      "This class holds several routines that access the data of multiple Atoms at once.")
      .def(py::init<std::vector<std::string>, Eigen::MatrixXd>())
      .def(py::init(&createGeometryPtr),
           "@brief The basic constructor using a python list of serenipy.Atom(s) or a .xyz file.\n"
           "Option 1: \n"
           "@param list of serenipy.Atom (internally std::shared_ptr<Atom>)\n"
           "Option2: \n"
           "@param string Path to an .xyz file.")
      .def("getCoordinates", &Geometry::getCoordinates,
           "@brief Getter for the coordinates.\n"
           "\n"
           "@returns MatrixXd containing all atom coordinates, dimensions: (nAtoms,3).")
      .def("getGradients", &Geometry::getGradients,
           "@brief Getter for the nuclear gradients if up-to-date.\n"
           "\n"
           "@returns MatrixXd containing all nuclear gradients, dimensions: (nAtoms,3).")
      .def("getNAtoms", &Geometry::getNAtoms, "@returns The number of atoms.")
      .def("getCoreCoreRepulsion", &Geometry::getCoreCoreRepulsion,
           "@returnsThe energy due to the Coulomb repulsion of the atoms from each other.")
      .def("getMaxX", &Geometry::getMaxX, "@returns the largest x coordinate of the underlying atoms")
      .def("getMaxY", &Geometry::getMaxY, "@returns the largest y coordinate of the underlying atoms")
      .def("getMaxZ", &Geometry::getMaxZ, "@returns the largest z coordinate of the underlying atoms")
      .def("getMinX", &Geometry::getMinX, "@returns the smallest (or most negative) x coordinate of the underlying atoms")
      .def("getMinY", &Geometry::getMinY, "@returns the smallest (or most negative) y coordinate of the underlying atoms")
      .def("getMinZ", &Geometry::getMinZ, "@returns the smallest (or most negative) z coordinate of the underlying atoms")
      .def("printGeometry", &Geometry::print, "@brief Prints the current geometry to the screen.")
      .def("printGradients", &Geometry::printGradients, "@brief Prints the current geometry gradients to the screen.")
      .def("setCoordinates", &Geometry::setCoordinates,
           "@brief Set all coordinates at once"
           "@param The new coordinates, dimensions: (nAtoms,3)")
      .def("setGradients", &Geometry::setGradients,
           "@brief Set the Gradients.\n"
           "@param The new geometry gradient, dimensions: (nAtoms,3)");
}