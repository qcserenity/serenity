/**
 * @file Triangle.cpp
 *
 * @date Feb 28, 2017
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

/* Include Class Header*/
#include "geometry/Triangle.h"
/* Include Serenity Internal Headers */
#include "geometry/Point.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

Triangle::Triangle(Point point1, Point point2, Point point3) : _sphereCenter(0, 0, 0), _radius(1.0) {
  _points.push_back(point1);
  _points.push_back(point2);
  _points.push_back(point3);
  _points.push_back(point1 + point2 + point3);
  _points[3] /= 3;
};

Triangle::Triangle(Point point1, Point point2, Point point3, Point center, double radius)
  : _sphereCenter(center), _radius(radius) {
  _points.push_back(point1);
  _points.push_back(point2);
  _points.push_back(point3);
  Point triangleCenter = (point1 + point2 + point3) / 3;
  Point move3 = triangleCenter - _sphereCenter;
  triangleCenter += move3 * (_radius / move3.distanceToOrigin() - 1.0);
  _points.push_back(triangleCenter);
}

double Triangle::getArea() const {
  auto ab = _points[1] - _points[0];
  auto ac = _points[2] - _points[0];

  double angle = ab.getX() * ac.getX() + ab.getY() * ac.getY() + ab.getZ() * ac.getZ();
  angle /= (ab.distanceToOrigin() * ac.distanceToOrigin());
  angle = acos(angle);

  double area = 0.5 * ab.distanceToOrigin() * ac.distanceToOrigin() * sin(angle);

  return area;
}

double Triangle::getArea(unsigned int calveLevel) const {
  double area = 0.0;
  if (calveLevel == 0) {
    area = this->getArea();
  }
  else {
    std::vector<Triangle> calvedTriangles = this->calve(true);
    for (auto tri : calvedTriangles)
      area += tri.getArea(calveLevel - 1);
  }
  return area;
}

const Eigen::Vector3d& Triangle::getNormVector() {
  if (!_normVector) {
    _normVector = std::make_shared<Eigen::Vector3d>(Eigen::Vector3d::Zero(3));
    double x = _sphereCenter.getX() - _points[3].getX();
    double y = _sphereCenter.getY() - _points[3].getZ();
    double z = _sphereCenter.getZ() - _points[3].getZ();
    *_normVector << x, y, z;
    *_normVector *= 1.0 / _normVector->norm();
  }
  return *_normVector;
}

std::vector<Triangle> Triangle::calve(bool fitToSurface) const {
  /*
   * get the new points P12=(P1+P2)/2
   */
  Point point12 = _points[0] + _points[1];
  Point point13 = _points[0] + _points[2];
  Point point23 = _points[1] + _points[2];
  point12 /= 2;
  point13 /= 2;
  point23 /= 2;

  if (fitToSurface) {
    Point move12 = point12 - _sphereCenter;
    Point move13 = point13 - _sphereCenter;
    Point move23 = point23 - _sphereCenter;
    point12 += move12 * (_radius / move12.distanceToOrigin() - 1.0);
    point13 += move13 * (_radius / move13.distanceToOrigin() - 1.0);
    point23 += move23 * (_radius / move23.distanceToOrigin() - 1.0);
  }

  /*
   * build up new triangles and put them in a vector
   */
  std::vector<Triangle> newTriangles = {Triangle(_points[0], point12, point13, _sphereCenter, _radius),
                                        Triangle(_points[1], point12, point23, _sphereCenter, _radius),
                                        Triangle(_points[2], point13, point23, _sphereCenter, _radius),
                                        Triangle(point12, point13, point23, _sphereCenter, _radius)};

  return newTriangles;
}

void Triangle::adjustToSphere(Point center, double radius) {
  _radius = radius;
  _sphereCenter = center;
  _points[3] /= (_points[3].distanceToOrigin());
  for (auto& point : _points) {
    point *= radius;
  }
  this->translate(center);
}

void Triangle::translate(Point translation) {
  for (auto& point : _points) {
    point += translation;
  }
}

} // namespace Serenity
