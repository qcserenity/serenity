/**
 * @file   Point.h
 * @author Thomas Dresselhaus
 *
 * @date   14. März 2014, 22:32
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
#ifndef POINT_H
#define POINT_H
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <array>
#include <cmath>

namespace Serenity {
/**
 * @class Point Point.h
 *
 * @brief x,y,z
 */
class Point : public std::array<double, 3> {
 public:
  /**
   * @param x
   * @param y
   * @param z
   */
  constexpr Point(double x, double y, double z) : std::array<double, 3>{x, y, z} {
  }

  constexpr Point(std::array<double, 3> point) : std::array<double, 3>(point) {
  }

  virtual ~Point() = default;

  virtual void addToX(double add_to_x) {
    this->at(0) += add_to_x;
  }

  virtual void addToY(double add_to_y) {
    this->at(1) += add_to_y;
  }

  virtual void addToZ(double add_to_z) {
    this->at(2) += add_to_z;
  }

  /// @returns the x-coordinate
  inline const double& getX() const {
    return this->at(0);
  }

  /// @param x new x coordinate
  virtual void setX(double x) {
    this->at(0) = x;
  }

  /// @returns the y-coordinate
  inline const double& getY() const {
    return this->at(1);
  }

  /// @param y new y coordinate
  virtual void setY(double y) {
    this->at(1) = y;
  }

  /// @returns the z-coordinate
  inline const double& getZ() const {
    return this->at(2);
  }

  /// @param z new z coordinate
  virtual void setZ(double z) {
    this->at(2) = z;
  }
  /**
   * @returns the absolute value of the vector (x,y,z)
   */
  double distanceToOrigin() const {
    return sqrt(this->at(0) * this->at(0) + this->at(1) * this->at(1) + this->at(2) * this->at(2));
  }
  /**
   * @param   rhs
   * @returns the shifted point
   */
  void operator+=(const Point& rhs) {
    this->at(0) += rhs[0];
    this->at(1) += rhs[1];
    this->at(2) += rhs[2];
  }
  /**
   * @param   rhs
   * @returns the shifted point
   */
  void operator-=(const Point& rhs) {
    this->at(0) -= rhs[0];
    this->at(1) -= rhs[1];
    this->at(2) -= rhs[2];
  }
  /**
   * @param   factor
   * @returns the shifted point
   */
  void operator*=(const double factor) {
    this->at(0) *= factor;
    this->at(1) *= factor;
    this->at(2) *= factor;
  }
  /**
   * @param   factor
   * @returns the shifted point
   */
  void operator/=(const double factor) {
    this->at(0) /= factor;
    this->at(1) /= factor;
    this->at(2) /= factor;
  }
  /**
   * @brief Operator ==
   * @param other The other point.
   * @return True, if the coordinates match.
   */
  bool operator==(const Point& other) const {
    return this->isSamePoint(other, 1e-9);
  }

  /**
   * @brief A function to check if two points are identical.
   * @param point The other point.
   * @param precision Precision to compare two doubles.
   * @return A bool determining if the two points are identical.
   */
  bool isSamePoint(const Point& point, double precision) const;

  /**
   * @brief Getter for the point as an Eigen::Vector3d.
   * @return
   */
  Eigen::Vector3d coords() const {
    return {this->getX(), this->getY(), this->getZ()};
  }
};

/**
 * @param   a, b
 * @returns the distance between a and b in atomic units
 */
inline double distance(const Point& a, const Point& b) {
  return sqrt((a.getX() - b.getX()) * (a.getX() - b.getX()) + (a.getY() - b.getY()) * (a.getY() - b.getY()) +
              (a.getZ() - b.getZ()) * (a.getZ() - b.getZ()));
}

Point operator+(Point lhs, const Point& rhs);
Point operator-(Point lhs, const Point& rhs);
Point operator*(Point lhs, const double rhs);
Point operator/(Point lhs, const double rhs);

} // namespace Serenity
#endif /* POINT_H */
