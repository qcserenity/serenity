/**
 * @file Plane.h
 *
 * @date   May 28, 2020
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

#ifndef GEOMETRY_PLANE_H_
#define GEOMETRY_PLANE_H_
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

class Line;

/**
 * @class Plane Plane.h
 * @brief A class that defines a plane in \f$ \mathbb{R}^3 \f$. A plane is given by all
 *        points \f$ \pmb{r} \in \mathbb{R}^3 \f$ that fulfill the equation
 *          \f$ 0 = (\pmb{r}-\pmb{p})\pmb{n}, \f$
 *        where \f$ \pmb{p} \f$ is a point in the plane, and  \f$ \pmb{n} \f$ is the vecotr
 *        orthogonal to the plane.
 */
class Plane {
 public:
  /**
   * @brief Constructor.
   * @param point The starting/reference point of the plane.
   * @param normalVector The vector orthogonal to the plane.
   */
  Plane(const Eigen::Vector3d& point, const Eigen::Vector3d normalVector);
  /**
   * @brief Getter for the reference point of the plane used for construction.
   */
  const Eigen::Vector3d& getPoint() const;
  /**
   * @brief Getter for the normal vector of the plane.
   */
  const Eigen::Vector3d& getNormalVector() const;
  /**
   * @brief Calculate the intersection between this plane and the plane other, with the
   *        constrain that the reference point of the resulting line is as close to the
   *        reference as possible.
   * @param other The other plane.
   * @param reference The reference point to which the starting point of the line will be
   *                  close to.
   * @param The intersection line between the planes. Returns a nullptr if the planes are
   *        parallel and not identical. If they are identical, a line with arbitrary direction
   *        and the starting point of this plane as a starting point. The constrain with respect
   *        to the reference is not applied!
   */
  std::shared_ptr<Line> calculateIntersection(const Plane& other, const Eigen::Vector3d& reference) const;
  /**
   * @brief Calculate the intersection between this plane and the plane other, with the
   *        constrain that the reference point of the resulting line is as close to the
   *        reference point of this planes as possible.
   * @param other The other plane.
   * @param The intersection line between the planes. Returns a nullptr if the planes are
   *        parallel and not identical. If they are identical, a line with arbitrary direction
   *        and the starting point of this plane as a starting point.
   */
  std::shared_ptr<Line> calculateIntersection(const Plane& other) const;
  /**
   * @brief Calculates the intersection point between a line and this plane.
   * @param line The line.
   * @return The intersection point if line and plane are intersecting, nullptr if not.
   *         If the line lies within the plane, the reference point of the line is returned.
   */
  std::shared_ptr<Eigen::Vector3d> calculateIntersection(const Line& line) const;

  /**
   * @brief Checks whether a given point is in the plane.
   * @param r The point.
   * @return True if in the plane. False otherwise.
   */
  bool pointIsInPlane(const Eigen::Vector3d r) const;
  /**
   * @brief Default destructor.
   */
  ~Plane();

 private:
  // The reference point.
  Eigen::Vector3d _p;
  // The normal vector.
  Eigen::Vector3d _n;
};

} /* namespace Serenity */

#endif /* GEOMETRY_PLANE_H_ */
