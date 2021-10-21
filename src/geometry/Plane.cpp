/**
 * @file Plane.cpp
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
#include "geometry/Plane.h"
/* Include Serenity Internal Headers */
#include "geometry/Line.h"

namespace Serenity {

Plane::Plane(const Eigen::Vector3d& point, const Eigen::Vector3d normalVector)
  : _p(point), _n(1.0 / normalVector.norm() * normalVector) {
}

const Eigen::Vector3d& Plane::getPoint() const {
  return _p;
}
const Eigen::Vector3d& Plane::getNormalVector() const {
  return _n;
}
std::shared_ptr<Line> Plane::calculateIntersection(const Plane& other) const {
  return calculateIntersection(other, _p);
}
std::shared_ptr<Line> Plane::calculateIntersection(const Plane& other, const Eigen::Vector3d& reference) const {
  const Eigen::Vector3d& nOther = other.getNormalVector();
  Eigen::Vector3d direction = _n.cross(nOther);
  double dirNorm = direction.norm();
  // Check if parallel.
  if (dirNorm < 1e-7) {
    // Check if the planes are identical.
    double testInPlane = (other.getPoint() - _p).dot(_n);
    if (std::abs(testInPlane) < 1e-7) {
      // Planes are identical. Generate any directional vector within the plane.
      const Eigen::Vector3d xAxis = {1, 0, 0};
      direction = _n.cross(xAxis);
      dirNorm = direction.norm();
      if (dirNorm < 1e-7) {
        const Eigen::Vector3d yAxis = {0, 1, 0};
        direction = _n.cross(yAxis);
        direction /= direction.norm();
      }
      return std::make_shared<Line>(_p, direction);
    }
    else {
      return nullptr;
    }
  }
  direction *= 1.0 / dirNorm;
  // Determine the line-point such that it is close to
  // the reference point and the points of both planes.
  // See:
  // https://math.stackexchange.com/questions/475953/how-to-calculate-the-intersection-of-two-planes
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(5, 5);
  M(0, 0) = 2.0;
  M(1, 1) = 2.0;
  M(2, 2) = 2.0;
  M.col(3).segment(0, 3) = _n;
  M.row(3).segment(0, 3) = _n;
  M.col(4).segment(0, 3) = nOther;
  M.row(4).segment(0, 3) = nOther;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(5);
  rhs << 2.0 * reference(0), 2.0 * reference(1), 2.0 * reference(2), _p.dot(_n), other.getPoint().dot(nOther);
  Eigen::VectorXd lagrange = M.inverse() * rhs;
  Eigen::Vector3d linePoint = lagrange.segment(0, 3);
  return std::make_shared<Line>(linePoint, direction);
}
std::shared_ptr<Eigen::Vector3d> Plane::calculateIntersection(const Line& line) const {
  const Eigen::Vector3d& d = line.getDirection();
  const Eigen::Vector3d& pl = line.getPoint();
  // Check if pl is in the plane
  if (std::abs((pl - _p).dot(_n)) < 1e-9)
    return std::make_shared<Eigen::Vector3d>(pl);
  // check if parallel.
  if (std::abs(d.dot(_n)) < 1e-9)
    return nullptr;
  // There must be a unique intersection.
  const double x = (_p - pl).dot(_n) / d.dot(_n);
  Eigen::Vector3d s = pl + x * d;
  return std::make_shared<Eigen::Vector3d>(s);
}

bool Plane::pointIsInPlane(const Eigen::Vector3d r) const {
  return std::abs((r - _p).dot(_n)) < 1e-9;
}

Plane::~Plane() = default;

} /* namespace Serenity */
