/**
 * @file MolecularSurface.h
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

#ifndef GEOMETRY_MOLECULARSURFACE_H_
#define GEOMETRY_MOLECULARSURFACE_H_

/* Include Serenity Internal Headers */
#include "grid/AtomCenteredGrid.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include <string>

namespace Serenity {

class MolecularSurfaceController;
class Sphere;

/**
 * @class
 * @brief A atom centered grid variant for a molecular surface.
 *
 * The surface is assumed to be constructed from a set of spheres.
 * To each sphere a set of points is associated.
 */
class MolecularSurface : public AtomCenteredGrid {
  friend MolecularSurfaceController;

 public:
  /**
   * @brief Constructor
   * @param gridPoints              The point coordinates.
   * @param weights                 The grid point weights.
   * @param normalVectors           The normal vectors.
   * @param label                   Characteristic label for the surface type.
   * @param sphereIndices           The point indices for each underlying sphere.
   * @param pointWiseSphereIndices  The sphere indices for each point.
   */
  MolecularSurface(std::unique_ptr<Eigen::Matrix3Xd> gridPoints, std::unique_ptr<Eigen::VectorXd> weights,
                   std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string label,
                   std::vector<std::pair<unsigned int, unsigned int>> sphereIndices,
                   std::vector<unsigned int> pointWiseSphereIndices, double r_s, std::vector<Sphere> spheres);
  /**
   * @brief Default destructor.
   */
  ~MolecularSurface();
  /**
   * @brief Getter for the normal vectors.
   * @return The normal vectors.
   */
  const Eigen::Matrix3Xd& getNormalVectors() const;
  /**
   * @brief Getter for the surface label.
   */
  std::string getLabel() const;
  /**
   * @brief Getter for the point wise sphere indices.
   * @return The point wise sphere indices.
   */
  const std::vector<unsigned int>& pointWiseSphereIndices() const;
  /**
   * @brief Getter for the probe radius used during the surface construction.
   * @return The probe radius.
   */
  double getRs() const;
  /**
   * @brief Getter for the spheres used during surface construction.
   * @return The spheres.
   */
  const std::vector<Sphere>& getSpheres() const;

 private:
  // The normal vectors.
  std::unique_ptr<Eigen::Matrix3Xd> _normalVectors;
  // The surface label.
  std::string _label;
  // The point to sphere mapping.
  std::vector<unsigned int> _pointWiseSphereIndices;
  // The probe radius.
  const double _r_s;
  // The spheres.
  const std::vector<Sphere> _spheres;
};

} /* namespace Serenity */

#endif /* GEOMETRY_MOLECULARSURFACE_H_ */
