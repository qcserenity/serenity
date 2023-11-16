/**
 * @file Line.cpp
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
/* Include Class Header*/
#include "geometry/Line.h"

namespace Serenity {

Line::Line(const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
  : _p(point), _d(1.0 / direction.norm() * direction) {
}

const Eigen::Vector3d& Line::getPoint() const {
  return _p;
}
const Eigen::Vector3d& Line::getDirection() const {
  return _d;
}
bool Line::isOnLine(const Eigen::Vector3d& r) const {
  const Eigen::Vector3d g = _p - r;
  double t = 0.0;
  if (std::abs(_d(0)) > 1e-9) {
    t = -g(0) / _d(0);
  }
  else if (std::abs(_d(1)) > 1e-9) {
    t = -g(1) / _d(1);
  }
  else if (std::abs(_d(2)) > 1e-9) {
    t = -g(2) / _d(2);
  }
  else {
    return g.norm() < 1e-9;
  }
  // Has to be the zero vector.
  Eigen::Vector3d test = g + t * _d;
  return test.norm() < 1e-9;
}

std::shared_ptr<Eigen::Vector3d> Line::getIntersection(const Line& other) const {
  const Eigen::Vector3d& p2 = other.getPoint();
  const Eigen::Vector3d& d2 = other.getDirection();
  if (this->isOnLine(p2))
    return std::make_shared<Eigen::Vector3d>(p2);
  const double f = 1.0 - (d2(0) * _d(1) / (_d(0) * d2(1)));
  const double t2 = 1.0 / f * ((_p(1) - p2(1)) / d2(1) + _d(1) * (p2(0) - _p(0)) / (_d(0) * d2(1)));
  const double t1 = (p2(0) - _p(0) + t2 * d2(0)) / _d(0);
  Eigen::Vector3d r = _p + t1 * _d;
  Eigen::Vector3d test = p2 + t2 * d2;
  const double diff = (r - test).norm();
  // Return nullptr if there is no intersection (points are not identical or nan elements).
  if (diff > 1e-9 || diff != diff)
    return nullptr;
  return std::make_shared<Eigen::Vector3d>(r);
}

Line::~Line() = default;

} /* namespace Serenity */
