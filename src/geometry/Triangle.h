/**
 * @file Triangle.h
 *
 * @date Feb 28, 2077
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

#ifndef GEOMETRY_TRIANGLE_H_
#define GEOMETRY_TRIANGLE_H_

/* Include Serenity Internal Headers */
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>
namespace Serenity {

class Point;

/*
 * A class resembling a triangle. This can be used to approximate
 * surfaces, enabling a discretization (e.g. for a COSMO algorithm)
 * or the calculation of surface areas.
 */
class Triangle {
 public:
  /**
   * @brief Constructor
   * @param point1-point3 The points of the triangle
   */
  Triangle(Point point1, Point point2, Point point3);

  Triangle(Point point1, Point point2, Point point3, Point center, double radius);

  virtual ~Triangle() = default;

  /**
   * @brief Getter function
   * @return The center of the triangle.
   */
  inline const Point& getCenter() const {
    return _points[3];
  }
  /**
   * @brief Getter for the triangle points. The fourth point is
   *        the triangle center.
   * @return The triangle points.
   */
  inline const std::vector<Point>& getPoints() const {
    return _points;
  }
  /**
   * @brief Calculates the area of one the triangle.
   * @return The area of the triangle.
   */
  double getArea() const;
  /**
   * @brief Getter for the area of the sphere covered by the triangle.
   * @param calveLevel The level to increase the resolution of the area calculation.
   * @return The covered area.
   */
  double getArea(unsigned int calveLevel) const;
  /**
   * @brief cuts the triangle in 4 smaller triangles by adding additional
   *        points in the center of each side.
   * @param fitToSurface If true, the triangles will be tilted to approximate
   *        a surface of a sphere with radius 1.
   * @return A vector containing the new triangles.
   */
  std::vector<Triangle> calve(bool fitToSurface) const;
  /**
   * @brief Adjusts the triangle to a sphere
   * @param center The center of the sphere.
   * @param radius The radius of the sphere.
   */
  void adjustToSphere(Point center, double radius);
  /**
   * @brief Scales the triangle about a factor
   * @param factor The factor to be scaled with.
   */
  void translate(Point translation);
  /**
   * @brief Getter for the normal vectors on the cavity.
   * @return
   */
  const Eigen::Vector3d& getNormVector();

 private:
  Point _sphereCenter;
  double _radius;
  std::vector<Point> _points;
  std::shared_ptr<Eigen::Vector3d> _normVector;
};

} // namespace Serenity

#endif /* BASICS_GEOMETRY_TRIANGLE_H_ */
