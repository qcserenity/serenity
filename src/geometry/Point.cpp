/**
 * @file   Point.cpp
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
/* Include Class Header*/
#include "geometry/Point.h"

namespace Serenity {
/**
 * @param lhs
 * @param rhs
 * @returns a point with coordinates which are the sum of the incoming points' coordinates
 */
Point operator+(Point lhs, const Point& rhs) {
  lhs += rhs;
  return lhs;
}
/**
 * @param lhs
 * @param rhs
 * @returns a point with coordinates which are the difference of the incoming points' coordinates
 */
Point operator-(Point lhs, const Point& rhs) {
  lhs -= rhs;
  return lhs;
}
/**
 * @param lhs
 * @param rhs
 * @returns a scaled version of the Point lhs
 */
Point operator*(Point lhs, const double rhs) {
  lhs *= rhs;
  return lhs;
}
/**
 * @param lhs
 * @param rhs
 * @returns a scaled version of the Point lhs
 */
Point operator/(Point lhs, const double rhs) {
  lhs /= rhs;
  return lhs;
}

bool Point::isSamePoint(const Point& point, double precision) const {
  double dX = std::abs(this->getX() - point.getX());
  double dY = std::abs(this->getY() - point.getY());
  double dZ = std::abs(this->getZ() - point.getZ());
  if (dX > precision || dY > precision || dZ > precision) {
    return false;
  }
  else {
    return true;
  }
}

} // namespace Serenity
