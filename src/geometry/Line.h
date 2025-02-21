/**
 * @file Line.h
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

#ifndef GEOMETRY_LINE_H_
#define GEOMETRY_LINE_H_
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/**
 * @class Line Line.h
 * @brief A class that defines a line. A line is given by the points \f$ \pmb(l)\f$ that fulfill:\n
 *    \f$ \pmb{l} = \pmb{p} + x \pmb{d}, x\in \mathbb{R}, \f$ \n
 *   where \f$ \pmb{p} \in \mathbb{R}^3 \f$ is some point on the line and \f$ \pmb{d}\in \mathbb{R}^3 \f$
 *   is the direction of the line.
 */
class Line {
 public:
  /**
   * @brief Constructor,
   * @param point The starting point.
   * @param direction The direction. The vector is normalized during construction.
   */
  Line(const Eigen::Vector3d& point, const Eigen::Vector3d& direction);
  /**
   * @brief Getter for the starting point vector.
   * @return The starting point vector.
   */
  const Eigen::Vector3d& getPoint() const;
  /**
   * @brief Getter for the direction vector.
   * @return The direction vector.
   */
  const Eigen::Vector3d& getDirection() const;
  /**
   * @brief Calculate the intersection between this line and another line.
   * @param other The other line.
   * @return The intersection point or nullptr if no intersection point exists.
   *		   If both lines are "identical" the starting point of other is returned.
   */
  std::shared_ptr<Eigen::Vector3d> getIntersection(const Line& other) const;
  /**
   * @brief Check if a point is on the line.
   * @param r The point.
   */
  bool isOnLine(const Eigen::Vector3d& r) const;
  /**
   * @brief Default destructor.
   */
  ~Line();

 private:
  // The reference/starting point.
  const Eigen::Vector3d _p;
  // The direction.
  const Eigen::Vector3d _d;
};

} /* namespace Serenity */

#endif /* GEOMETRY_LINE_H_ */
