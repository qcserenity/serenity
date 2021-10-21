/**
 * @file MolecularSurface.cpp
 *
 * @date   Nov 12, 2020
 * @author Moritz Bensberg
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

/* Include Class Header*/
#include "geometry/MolecularSurface.h"
/* Include Serenity Internal Headers */
#include "geometry/Sphere.h"

namespace Serenity {

MolecularSurface::MolecularSurface(std::unique_ptr<Eigen::Matrix3Xd> gridPoints, std::unique_ptr<Eigen::VectorXd> weights,
                                   std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string label,
                                   std::vector<std::pair<unsigned int, unsigned int>> sphereIndices,
                                   std::vector<unsigned int> pointWiseSphereIndices, double r_s, std::vector<Sphere> spheres)
  : AtomCenteredGrid(std::move(gridPoints), std::move(weights), sphereIndices),
    _normalVectors(std::move(normalVectors)),
    _label(label),
    _pointWiseSphereIndices(pointWiseSphereIndices),
    _r_s(r_s),
    _spheres(spheres) {
  this->_sorted = false;
}

const Eigen::Matrix3Xd& MolecularSurface::getNormalVectors() const {
  return *_normalVectors;
}

std::string MolecularSurface::getLabel() const {
  return _label;
}

const std::vector<unsigned int>& MolecularSurface::pointWiseSphereIndices() const {
  return _pointWiseSphereIndices;
}

double MolecularSurface::getRs() const {
  return _r_s;
}

const std::vector<Sphere>& MolecularSurface::getSpheres() const {
  return _spheres;
}

MolecularSurface::~MolecularSurface() = default;
} /* namespace Serenity */
