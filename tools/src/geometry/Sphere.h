/**
 * @file Sphere.h
 *
 * @date Sep 23, 2016
 * @author David Schnieders
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

#ifndef GEOMETRY_SPHERE_H_
#define GEOMETRY_SPHERE_H_
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
#include "geometry/Triangle.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class
 * @brief The different sphere types.
 *   VanDerWaals: Initial sphere constructed from an atom.
 *   Ghost: Sphere that is not centered on an atom.
 *   ZeroArea: Sphere thats surface is fully covered by other spheres.
 *   Engulfed: Sphere that is completely within other spheres.
 */
enum class SphereType { VanDerWaals, Ghost, ZeroArea, Engulfed };
/**
 * @class
 * @brief A class that represents a sphere for molecular cavity construction.
 */
class Sphere {
 public:
  /**
   * @brief Constructor
   * @param center The center of the sphere.
   * @param radius The radius of the sphere.
   * @param sphereType The sphere type.
   * @param angularMomentum The angular momentum.
   */
  Sphere(Point center, double radius, SphereType sphereType, unsigned int angularMomentum);
  virtual ~Sphere() = default;
  /**
   * @brief Getter for the sphere center.
   * @return The sphere center.
   */
  Point getCenter() const {
    return _center;
  }
  /**
   * @brief Getter for the sphere center.
   * @return The sphere center as an Eigen::Vector3d.
   */
  inline const Eigen::Vector3d& getCenterCoords() const {
    return _centerCoords;
  }
  /**
   * @brief Getter for the sphere radius.
   * @return The sphere radius.
   */
  inline const double& getRadius() const {
    return _radius;
  }
  /**
   * @brief Getter for the sphere type.
   * @return The sphere type.
   */
  SphereType getSphereType() const {
    return _sphereType;
  }
  /**
   * @brief Setter for the sphere center.
   * @param newCenter The new sphere center.
   */
  void setCenter(Point newCenter) {
    _center = newCenter;
    _centerCoords << _center.getX(), _center.getY(), _center.getZ();
  }
  /**
   * @brief Setter for the sphere radius.
   * @param newRadius The new radius.
   */
  void setRadius(double newRadius) {
    _radius = newRadius;
  }
  /**
   * @brief Setter for the sphere type.
   * @param newSphereType The new sphere type.
   */
  void setSphereType(SphereType newSphereType) {
    _sphereType = newSphereType;
  }
  /**
   * @brief Getter for angular momentum used for the Delley-type cavity.
   * @return The angular momentum.
   */
  unsigned int getAngularMomentum() const {
    return _angularMomentum;
  }
  /**
   * @brief Calculates the distance between the centers of two spheres
   * @param otherSphere The other sphere the distance should be calculated to
   * @return The distance between the centers of the two spheres.
   */
  double distanceTo(Sphere otherSphere) const {
    double distance = (_center - otherSphere.getCenter()).distanceToOrigin();
    return distance;
  }

  /**
   * @brief The solvent accessible surface of the sphere is approximated
   *        with 60 triangles. A sphere of radius 1 is approximated first
   *        and saved as static object to prevent multiple calculations
   *        of this step. The resulting triangles are then scaled to the
   *        correct radius (i.e. sphereRadius+solventRadius) and returned.
   * @param solventRadius The radius of the solvent.
   * @return A vector containing the triangles which approximate the surface.
   */
  std::vector<Triangle> getSolventAccessibleSurface(double solventRadius);

  /**
   * @brief Approximates the surface of the sphere with 60
   *        triangles
   * @return A vector containing the triangles which approximate the surface.
   */
  std::vector<Triangle> getSurface() {
    auto surface = getSolventAccessibleSurface(0.0);
    return surface;
  }

 private:
  /**
   * @brief Builds up 60 triangles to approximate the surface of the Sphere.
   *
   * Adopted from the GEPOL algorithm:
   * http://www.ccl.net/cca/software/SOURCES/FORTRAN/molecular_surface/gepol93\n\n
   *
   * Ref: J Comp. Chem. 8, 778--787 (1987)
   */
  void calcSurface();
  // Center of the sphere.
  Point _center;
  // Center of the sphere as an reasonable vector.
  Eigen::Vector3d _centerCoords;
  // Radius of the sphere.
  double _radius;
  // The sphere type.
  SphereType _sphereType;
  // The angular momentum.
  unsigned int _angularMomentum;
  /*
   * A static container for the surface of a sphere with radius 1
   * This only needs to be calculated once and can then be used by
   * every Sphere object
   */
  static std::vector<Triangle> _normSurface;
};

} // namespace Serenity

#endif /* BASICS_GEOMETRY_SPHERE_H_ */
