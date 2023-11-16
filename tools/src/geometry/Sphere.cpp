/**
 * @file Sphere.cpp
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
#include "geometry/Sphere.h"

namespace Serenity {

std::vector<Triangle> Sphere::_normSurface = {};

Sphere::Sphere(Point center, double radius, SphereType sphereType, unsigned int angularMomentum)
  : _center(center), _radius(radius), _sphereType(sphereType), _angularMomentum(angularMomentum) {
  _centerCoords << _center.getX(), _center.getY(), _center.getZ();
};

std::vector<Triangle> Sphere::getSolventAccessibleSurface(double solventRadius) {
  // a calculated _normSurface should be of length 60
  if (_normSurface.size() != 60) {
    calcSurface();
  }
  auto surface = _normSurface;
  for (auto& triangle : surface) {
    triangle.adjustToSphere(_center, _radius + solventRadius);
  }

  return surface;
}

void Sphere::calcSurface() {
  /*
   * The points which form each triangle of the surface.
   */
  std::vector<std::vector<unsigned int>> triangleMap = {
      {0, 5, 1},    {0, 1, 2},    {0, 2, 3},    {0, 3, 4},    {0, 4, 5},    {6, 1, 5},    {7, 2, 1},    {8, 3, 2},
      {9, 4, 3},    {10, 5, 4},   {7, 1, 11},   {8, 2, 12},   {9, 3, 13},   {10, 4, 14},  {6, 5, 15},   {6, 11, 1},
      {7, 12, 2},   {8, 13, 3},   {9, 14, 4},   {10, 15, 5},  {7, 11, 17},  {8, 12, 18},  {9, 13, 19},  {10, 14, 20},
      {6, 15, 16},  {6, 16, 11},  {7, 17, 12},  {8, 18, 13},  {9, 19, 14},  {10, 20, 15}, {21, 11, 16}, {22, 12, 17},
      {23, 13, 18}, {24, 14, 19}, {25, 15, 20}, {21, 17, 11}, {22, 18, 12}, {23, 19, 13}, {24, 20, 14}, {25, 16, 15},
      {21, 16, 26}, {22, 17, 27}, {23, 18, 28}, {24, 19, 29}, {25, 20, 30}, {21, 27, 17}, {22, 28, 18}, {23, 29, 19},
      {24, 30, 20}, {25, 26, 16}, {21, 26, 27}, {22, 27, 28}, {23, 28, 29}, {24, 29, 30}, {25, 30, 26}, {31, 27, 26},
      {31, 28, 27}, {31, 29, 28}, {31, 30, 29}, {31, 26, 30}};
  /*
   * The angles used to calculate the coordinates of the points on the surface.
   * The underlying pentakisdodecahedron can be constructed from 6 pentagons
   * above one another plus one point above the top and one point below the
   * lowest pentagon.
   *  Schematic side view:
   *                          *             top
   *                       *     *          1.
   *                    *          *        2.
   *                  *              *      3.
   *                    *          *        4.
   *                       *     *          5.
   *                          *             bottom.
   * Here 'thev' denotes the angle for deviation from the z-axis.
   *
   * The pentagons have either their first point on the x-axis or
   * the pentagon is rotated by pi/5. This off-set is contained in
   * 'fiv'.
   * The points of each pentagon lie in the same xy-plane where the
   * angle between each point is 2pi/5:
   *
   *                     *   *
   *          2pi/5 --->  \
   *                  *----x    *
   *
   *                       *
   * This is contained in fir.
   * All points are now constructed by looping over the pentagons, applying the
   * starting offset for the deviation from the y-axis and incrementing the angle
   * with 2pi/5.
   */
  std::vector<double> thev = {0.6523581397843682, 1.1071487177940905, 1.3820857960113345,
                              1.7595068575784587, 2.0344439357957027, 2.4892345138054251};
  std::vector<double> fiv = {0.6283185307179586, 0.0, 0.6283185307179586, 0.0, 0.6283185307179586, 0.0};
  double fir = 1.2566370614359173;

  /*
   * Build up the points which will be used to generate the triangles later
   * They are calculated by simple trigonometric functions, the corresponding
   * angles can be found hard coded above.
   */
  std::vector<Point> trianglePoints;
  trianglePoints.push_back(Point(0.0, 0.0, 1.0));
  for (unsigned int i = 0; i < 6; i++) {
    double cosTh = cos(thev[i]);
    double sinTh = sin(thev[i]);
    double fivTmp = fiv[i];
    for (unsigned int j = 0; j < 5; j++) {
      fivTmp += fir;
      if (j == 0) {
        fivTmp = fiv[i];
      }
      Point tmp(sinTh * cos(fivTmp), sinTh * sin(fivTmp), cosTh);
      trianglePoints.push_back(tmp);
    }
  }
  trianglePoints.push_back(Point(0.0, 0.0, -1.0));

  /*
   * build up the triangles using the points generated above.
   */
  for (unsigned int triangle = 0; triangle < triangleMap.size(); triangle++) {
    Triangle tmpTriangle(trianglePoints[triangleMap[triangle][0]], trianglePoints[triangleMap[triangle][1]],
                         trianglePoints[triangleMap[triangle][2]]);
    _normSurface.push_back(tmpTriangle);
  }
}

} // namespace Serenity
