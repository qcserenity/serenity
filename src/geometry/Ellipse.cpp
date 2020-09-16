/**
 * @file Ellipse.cpp
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
#include "geometry/Ellipse.h"
/* Include Serenity Internal Headers */
#include "geometry/Line.h"
#include "geometry/Plane.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

Ellipse::Ellipse(const Eigen::Vector3d& center, const Eigen::Vector3d& rad1, const Eigen::Vector3d& rad2)
  : _center(center), _rad1(rad1), _rad2(rad2) {
  if (std::abs(_rad1.dot(_rad2)) > 1e-9)
    throw SerenityError("The radii vectors of the ellipse have to be orthogonal");
}

double Ellipse::getArea() {
  return M_PI * _rad1.norm() * _rad2.norm();
}

const Eigen::Vector3d& Ellipse::getR1() {
  return _rad1;
}
const Eigen::Vector3d& Ellipse::getR2() {
  return _rad2;
}
const Eigen::Vector3d& Ellipse::getCenter() {
  return _center;
}
const Plane& Ellipse::getPlane() {
  if (!_plane) {
    Eigen::Vector3d n = _rad1.cross(_rad2);
    n *= 1.0 / n.norm();
    _plane = std::make_shared<Plane>(_center, n);
  }
  return *_plane;
}
void Ellipse::scaleByArea(double newArea) {
  const double radiiScaling = sqrt(newArea / this->getArea());
  _rad1 *= radiiScaling;
  _rad2 *= radiiScaling;
}
void Ellipse::scaleByFactor(double factor) {
  const double radiiScaling = sqrt(factor);
  _rad1 *= radiiScaling;
  _rad2 *= radiiScaling;
}
std::vector<Eigen::Vector3d> Ellipse::getIntersectionPoints(const Line& line) {
  std::vector<Eigen::Vector3d> intersections = {};
  const Eigen::Vector3d& pl = line.getPoint();
  const Eigen::Vector3d& pe = _center;
  const Eigen::Vector3d& d = line.getDirection();
  unsigned int ui = 1; // lower and upper index in devision below.
  unsigned int li = 0; // One entry of d has to be != 0. Otherwise the line
                       // construction should have failed.
  if (std::abs(d(li)) < 1e-6) {
    li = 2;
    if (std::abs(d(li)) < 1e-6) {
      li = 1;
      ui = 0;
    }
  }
  const double dydx = d(ui) / d(li);
  const double c = pl(ui) - pe(ui) + dydx * (pe(li) - pl(li));
  const double fc = _rad1(ui) - dydx * _rad1(li);
  const double fs = _rad2(ui) - dydx * _rad2(li);
  const double argument = fc * fc + fs * fs - c * c;
  const double fcPlusc = fc + c;
  if (std::abs(fcPlusc) > 1e-6 && argument > 1e-9) {
    const double t1 = 2.0 * atan((fs + sqrt(argument)) / fcPlusc);
    const double t2 = 2.0 * atan((fs - sqrt(argument)) / fcPlusc);
    const Eigen::Vector3d s1 = pe + cos(t1) * _rad1 + sin(t1) * _rad2;
    const Eigen::Vector3d s2 = pe + cos(t2) * _rad1 + sin(t2) * _rad2;
    intersections.push_back(s1);
    intersections.push_back(s2);
  }
  else if (std::abs(fcPlusc) <= 1e-6 && std::abs(fs) > 1e-6) {
    const double t1 = -2.0 * atan(fc / fs);
    const double t2 = M_PI;
    const Eigen::Vector3d s1 = pe + cos(t1) * _rad1 + sin(t1) * _rad2;
    intersections.push_back(s1);
    if (std::abs(t1 - t2) > 1e-6) {
      const Eigen::Vector3d s2 = pe + cos(t2) * _rad1 + sin(t2) * _rad2;
      intersections.push_back(s2);
    }
  }
  return intersections;
}
bool Ellipse::cutWithPlane(const Plane& plane, const Eigen::Vector3d& reference, CutEllipse& rFragment, CutEllipse& dFragment) {
  bool cutsPlane = true;
  double lArea = this->getArea();
  double rArea = 0.0;
  Eigen::Vector3d rCenter = Eigen::Vector3d::Zero(3);
  Eigen::Vector3d lCenter = _center;
  // Check the normal vectors of the plane to have an early bail out for parallel/identical planes.
  const Eigen::Vector3d& ellipsePlaneNormal = this->getPlane().getNormalVector();
  bool isParallel = std::abs(ellipsePlaneNormal.cross(plane.getNormalVector()).norm()) < 1e-9;

  // Get/check intersection between planes.
  std::shared_ptr<Line> intersection = nullptr;
  if (!isParallel)
    intersection = this->getPlane().calculateIntersection(plane, _center);
  if (intersection) {
    // Get intersection between ellipse and line.
    std::vector<Eigen::Vector3d> intersectionPoints = this->getIntersectionPoints(*intersection);
    // if the line and the ellipse do not touch or the line is a tangent.
    cutsPlane = intersectionPoints.size() > 1;
    if (cutsPlane) {
      // If reached here, the ellipse is cut by the plane.
      // Transform the ellipse to the origin and rotate it such that
      // the _rad1 is on the x-axis and _rad2 is the y-axis.
      double a = _rad1.norm();
      double b = _rad2.norm();
      Eigen::Vector3d xAxis;
      Eigen::Vector3d yAxis;
      xAxis << 1, 0, 0;
      yAxis << 0, 1, 0;
      // 1. Rotate around normal vector on _rad1 and the x-axis.
      Eigen::Vector3d nR1X = _rad1.cross(xAxis);
      const double nR1XNorm = nR1X.norm();
      Eigen::Matrix3d toXAxis = Eigen::Matrix3d::Identity(3, 3);
      Eigen::Matrix3d toXAxis_inv = Eigen::Matrix3d::Identity(3, 3);
      if (nR1XNorm > 1e-6) {
        nR1X *= 1.0 / nR1XNorm;
        const double angleToX = acos(_rad1(0) / a); // Dot-product with x-axis is trivial.
        toXAxis = getRotationMatrixinR3(nR1X, angleToX);
        toXAxis_inv = getRotationMatrixinR3(nR1X, -angleToX);
      }
      // Rotate _rad2.
      const Eigen::Vector3d rad2Rot = toXAxis * _rad2;
      // 2. Rotate the "new rad2" to the y-axis by rotating around the x-axis.
      double angleToY = acos(rad2Rot(1) / b);
      // Ensure right handed coordinate system for the rotation.
      Eigen::Vector3d rotationAxis = rad2Rot.cross(yAxis);
      double rotationAxisNorm = rotationAxis.norm();
      Eigen::Matrix3d toYAxis = Eigen::MatrixXd::Identity(3, 3);
      Eigen::Matrix3d toYAxis_inv = Eigen::MatrixXd::Identity(3, 3);
      if (rotationAxisNorm > 1e-6) {
        rotationAxis /= rotationAxis.norm();
        toYAxis = getRotationMatrixinR3(rotationAxis, angleToY);
        toYAxis_inv = getRotationMatrixinR3(rotationAxis, -angleToY);
      }
      // 3. Rotate the intersection points.
      for (auto& point : intersectionPoints) {
        point = toYAxis * toXAxis * (point - _center);
        //				if(std::abs(point(2)) > 1e-6 ) throw SerenityError((std::string)"This point has
        // to be in the xy-plane. Z-coordinate "+point(2));
      }
      Eigen::Vector3d C = intersectionPoints[1];
      Eigen::Vector3d B = intersectionPoints[0];
      // Make sure that xB >= xC. If not. Switch B and C.
      if (B(0) < C(0)) {
        const Eigen::Vector3d tmp = B;
        B = C;
        C = tmp;
      }
      double deltaX = B(0) - C(0);
      double deltaY = B(1) - C(1);
      double xA = 0.0;
      if (std::abs(deltaX) > 1e-6) {
        // Check if y-coordinates are identical or the cut is below the x-axis:
        /*
         *                          |
         *                      *   |b  *
         *              *           |            *
         *         *        uEArea  |                 *
         *   -a *                   |                    * a
         *  --*---------------------O--------A-------------*-->
         *      *                   |                    *
         *         *        lEArea  |                 *
         *       -------C-------------------------B------
         *                       *  |-b  *
         *                          |
         * If such a case occurs: Switch x and y axis by rotation by pi/2.
         *
         * Check if there is a negative slope between B and C.. If so, invert the x-axis.
         *
         *                          |
         *                      *   |b  B
         *              *           |    \       *
         *         *        uEArea  |     \           *
         *   -a *                   |      \             * a
         *  --*---------------------O-------A--------------*-->
         *      *                   |        \           *
         *         *        lEArea  |         \       *
         *              *           |          C  *
         *                       *  |-b  C
         *                          |
         *
         * Otherwise, the assumptions below will not hold.
         */

        if (std::abs(deltaY) < 1e-6 || (C(1) < 0 && B(1) < 0) || (C(1) > 0 && B(1) > 0)) {
          const Eigen::MatrixXd rotatePiHalf = getRotationMatrixinR3({0, 0, 1}, M_PI / 2.0);
          const Eigen::MatrixXd rotatePiHalf_inv = getRotationMatrixinR3({0, 0, 1}, -M_PI / 2.0);
          toYAxis = rotatePiHalf * toYAxis;
          toYAxis_inv = toYAxis_inv * rotatePiHalf_inv;
          C = rotatePiHalf * C;
          B = rotatePiHalf * B;
          // switch a and b.
          const double tmp = a;
          a = b;
          b = tmp;
        }
        deltaY = B(1) - C(1);
        deltaX = B(0) - C(0);
        double slopeCB = deltaY / deltaX;
        if (slopeCB < 0) {
          Eigen::MatrixXd invertX(3, 3);
          invertX << -1, 0, 0, 0, 1, 0, 0, 0, 1;
          toYAxis = invertX * toYAxis;
          toYAxis_inv = toYAxis_inv * invertX;
          C = invertX * C;
          B = invertX * B;
          slopeCB *= -1.0;
        }
        if (std::abs(deltaX) > 1e-6) {
          xA = C(0) - C(1) / slopeCB; // yC + s * dx = 0 <=>dx = -yC/s => x = xC+dx.
        }
        else {
          xA = B(0);
        }
      }
      else {
        xA = B(0);
        if (B(1) < C(1)) {
          const Eigen::Vector3d tmp = B;
          B = C;
          C = tmp;
        }
      }

      // Make sure that xB >= xC. If not. Switch B and C.
      if (B(0) < C(0)) {
        const Eigen::Vector3d tmp = B;
        B = C;
        C = tmp;
      }
      /*
       * We can calculate the area of both fragments by splitting
       * one fragment into two triangles and two functions under the
       * ellipse:
       *                          |
       *                      *   |b  *        /
       *              *	          |    	      B *
       *         *        uEArea  |          /|     *
       *   -a *                   |         / |        * a
       *  --*---------------------O--------A-------------*-->
       *      *                   |    |  /            *
       *         *        lEArea  |    | /          *
       *              *           |    |/       *
       *                       *  |-b  C
       *	                        |   /
       * The area of the left fragment is then obtained by summing
       * over the fragment areas and subtracting the triangle that contains
       * A and B.
       * The area of the second fragment is then obtained from the total area
       * minus the left fragment's area.
       */
      bool xABetween_xB_xC = C(0) <= xA && xA <= B(0);
      // Ignore cases where the ellipse is only just cut in a tiny part.
      if (xABetween_xB_xC) {
        const double xB = B(0);
        const double yB = B(1);
        const double xC = C(0);
        const double yC = C(1);
        const double uEArea = integrate(xB, a, b);
        const double lEArea = integrate(xC, a, b);
        const double uTriArea = 0.5 * (xB - xA) * std::abs(yB);
        const double lTriArea = 0.5 * (xA - xC) * std::abs(yC);
        const double leftArea = uEArea + lEArea - uTriArea + lTriArea;
        const double rightArea = lArea - leftArea;
        /*
         * The center of gravity is given by G = 1/A * int_F r dr,
         * where the integral runs over the area of function F.
         * Since we have transformed the ellipse to the origin, we now that
         * int_F r dr = 0 for F = E (E being the ellipse) (*).
         * It is sufficient to calculate the integral over one fragment, use (*)
         * to obtain the integral for the other fragment, and then weight by the area.
         * As before, we split our fragment into two areas below the ellipse
         * and two triangles.
         */
        if (leftArea > 1e-9 && rightArea > 1e-9) {
          Eigen::Vector3d uTriCGrav = B;
          uTriCGrav(0) += xA + xB;
          uTriCGrav *= uTriArea / 3.0;
          Eigen::Vector3d lTriCGrav = C;
          lTriCGrav(0) += xA + xC;
          lTriCGrav *= lTriArea / 3.0;
          const Eigen::Vector3d uECGrav = centerOfGravityEllipseFragment(xB, a, b, +1);
          const Eigen::Vector3d lECGrav = centerOfGravityEllipseFragment(xC, a, b, -1); // Flip y-axis. Thus, the sign
                                                                                        // is -1.
          Eigen::Vector3d leftCenterGrav = uECGrav + lECGrav - uTriCGrav + lTriCGrav;
          Eigen::Vector3d rightCenterGrav = -1.0 / rightArea * leftCenterGrav;
          leftCenterGrav = 1.0 / leftArea * leftCenterGrav;
          // Transform back into the original frame.
          leftCenterGrav = toXAxis_inv * toYAxis_inv * leftCenterGrav + _center;
          rightCenterGrav = toXAxis_inv * toYAxis_inv * rightCenterGrav + _center;
          // Save results.
          lArea = leftArea;
          lCenter = leftCenterGrav;
          rArea = rightArea;
          rCenter = rightCenterGrav;
        }
        else if (leftArea > 1e-9) {
          lArea = this->getArea();
          lCenter = this->getCenter();
          rArea = 0.0;
          rCenter = Eigen::Vector3d::Zero(3);
        }
        else {
          rArea = this->getArea();
          rCenter = this->getCenter();
          lArea = 0.0;
          lCenter = Eigen::Vector3d::Zero(3);
        } // else leftArea > 1e-9 && rightArea > 1e-9
      }   // if xABetween_xB_xC
      else {
        // Tiny cut. Assume no cut!
        cutsPlane = false;
      }
    }
  } // if !isParallel
  // Check if the planes are identical
  if (plane.pointIsInPlane(_center) && isParallel) {
    rArea = 0.0;
    lArea = 0.0;
    rCenter.setZero();
    lCenter.setZero();
  }
  // Check on which side of the plane the fragments are.
  // In order to do so, we need to check if lCenter and the reference
  // are 'above' or below the plane. This coordinates axis, orthogonal
  // to the plane, happens to be the normal vector of the plane.
  // We then check the sign of the coordinate with respect to its orthogonal
  // projection to the plane:
  /*
   *      <--------------- normal vector ----
   *           l--->|
   *                |
   *                |<--o
   *      r-------->|
   *                |
   *            the plane
   */
  const Line rayToPlaneL(lCenter, plane.getNormalVector());
  std::shared_ptr<Eigen::Vector3d> sL = plane.calculateIntersection(rayToPlaneL);
  const Line rayToPlaneRef(reference, plane.getNormalVector());
  std::shared_ptr<Eigen::Vector3d> sR = plane.calculateIntersection(rayToPlaneRef);

  const double normalCoordRef = (reference - *sR).dot(plane.getNormalVector());
  const double normalCoordL = (lCenter - *sL).dot(plane.getNormalVector());

  bool onReferenceSide = true;
  if ((normalCoordRef >= 0 && normalCoordL < 0) || (normalCoordRef < 0 && normalCoordL >= 0))
    onReferenceSide = false;
  if (not onReferenceSide) {
    double tmpArea = rArea;
    rArea = lArea;
    lArea = tmpArea;
    const Eigen::Vector3d tmpCenter = rCenter;
    rCenter = lCenter;
    lCenter = tmpCenter;
  }
  // Assign to final objects and return.
  rFragment.area = lArea;
  rFragment.centerOfGravity = lCenter;
  dFragment.area = rArea;
  dFragment.centerOfGravity = rCenter;
  return not isParallel && cutsPlane;
}

double Ellipse::integrate(double x, double a, double b) {
  const double sqroot = sqrt(a * a - x * x);
  const double func_Ma = -b * a * M_PI / 4.0;
  return b * (x * sqroot + a * a * atan(x / sqroot)) / (2.0 * a) - func_Ma;
}
double Ellipse::integrateXY(double x, double a, double b) {
  const double aaxx = a * a - x * x;
  const double cubed = aaxx * aaxx * aaxx;
  return -b * sqrt(cubed) / (3.0 * a);
}
double Ellipse::integrateYY(double x, double a, double b) {
  const double func_Ma = -2.0 / 3.0 * a * b * b;
  return 0.5 * (b * b * (x - x * x * x / (3.0 * a * a)) - func_Ma);
}
Eigen::Vector3d Ellipse::centerOfGravityEllipseFragment(double end, double a, double b, int ySign) {
  const double r_x = integrateXY(end, a, b);
  const double r_y = ySign * integrateYY(end, a, b);
  Eigen::Vector3d r;
  r << r_x, r_y, 0.0;
  return r;
}

Eigen::Matrix3d Ellipse::getRotationMatrixinR3(Eigen::Vector3d n, double a) {
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero(3, 3);
  const double cos_a = cos(a);
  const double sin_a = sin(a);
  const double e_cos_a = 1.0 - cos_a;
  // clang-format off
  R << n(0)*n(0) * e_cos_a+cos_a,         n(0) * n(1) * e_cos_a - n(2)*sin_a,       n(0) * n(2) * e_cos_a + n(1) * sin_a,
       n(1)*n(0) * e_cos_a+n(2) * sin_a,  n(1) * n(1) * e_cos_a + cos_a,            n(1) * n(2) * e_cos_a - n(0) * sin_a,
       n(2)*n(0) * e_cos_a-n(1) * sin_a,  n(2) * n(1) * e_cos_a + n(0)*sin_a,       n(2) * n(2) * e_cos_a + cos_a;
  // clang-format on
  return R;
}

Ellipse::~Ellipse() = default;

} /* namespace Serenity */
