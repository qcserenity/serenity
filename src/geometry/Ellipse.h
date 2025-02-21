/**
 * @file Ellipse.h
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

#ifndef GEOMETRY_ELLIPSE_H_
#define GEOMETRY_ELLIPSE_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

class Plane;
class Line;
/**
 * @struct
 * @brief A struct that defines the fragments of an ellipse.
 */
struct CutEllipse {
  double area = 0.0;
  Eigen::Vector3d centerOfGravity = Eigen::Vector3d::Zero(3);
};

/**
 * @class Ellipse Ellipse.h
 * @brief A class that models an ellipse.
 *   The ellipse is defined by the points \f$ \pmb{r} \in \mathbb{R}^3 \f$ that fulfill the equation\n
 *     \f$ \pmb{r} = \pmb{c} + \cos(t) \pmb{r}_1 + \sin(t) \pmb{r}_2 \f$\n
 *   for any \f$ t \in \mathbb{R}\f$. The center of the ellipse is denoted by \f$ \pmb{c} \f$ and
 *   \f$ \pmb{r}_1 \f$ and \f$ \pmb{r}_2 \f$ are translation vectors.
 */
class Ellipse {
 public:
  /**
   * @brief Constructor.
   * @param center The center of the ellipse.
   * @param rad1 The first translation vector.
   * @param rad2 The second translation vector (has to be orthogonal to rad1).
   */
  Ellipse(const Eigen::Vector3d& center, const Eigen::Vector3d& rad1, const Eigen::Vector3d& rad2);
  /**
   * @brief Getter for the area of the ellipse.
   * @return The area.
   */
  double getArea();
  /**
   * @brief Getter for the translation vector r1.
   * @return The translation vector r1.
   */
  const Eigen::Vector3d& getR1();
  /**
   * @brief Getter for the translation vector r2.
   * @return The translation vector r2.
   */
  const Eigen::Vector3d& getR2();
  /**
   * @brief Getter for the center of the ellipse.
   * @return
   */
  const Eigen::Vector3d& getCenter();
  /**
   * @brief Getter for the plane in which the ellipse is in.
   * @return The plane.
   */
  const Plane& getPlane();
  /**
   * @brief Calculate the intersection points of the line with the ellipse.
   * @param line The intersecting line.
   * @return The intersection points. If the line is touching or passing the ellipse
   *         the vector will have the length of 0. Otherwise it will contain the pair
   *         of intersecting points.
   *         Accuracy of evaluating whether it is passing should be around 1e-9.
   */
  std::vector<Eigen::Vector3d> getIntersectionPoints(const Line& line);
  /**
   * @brief Cut the ellipse with a plane and calculate the center of gravity of the fragments and their areas.
   *        The areas and centers are updated within the fragments parsed to this function.
   *        If the cutting plane is identical with the plane of the ellipse, both areas and centers are set to
   *        zero.
   *        The rFragment will always describe the fragment close to the reference point and the dFragment the
   *        remaining part. One of the fragments may be set to zero if the plane fully hides the ellipse from the
   *        reference or the ellipse is on the side of the reference and the plane does not cut it.
   * @param plane     The cutting plane.
   * @param reference The reference point.
   * @param rFragment The fragment close to the reference point.
   * @param dFragment The fragment further away from the reference point.
   * @return True if the plane cuts the ellipse and the planes are not identical. False otherwise.
   */
  bool cutWithPlane(const Plane& plane, const Eigen::Vector3d& reference, CutEllipse& rFragment, CutEllipse& dFragment);
  /**
   * @brief Scale the radii of the ellipse such that the ellipse will have the given area.
   * @param newArea The new area.
   */
  void scaleByArea(double newArea);
  /**
   * @brief Scale the radii of the ellipse such that the total area of it is scaled by the given factor.
   * @param factor The scaling factor.
   */
  void scaleByFactor(double factor);
  /**
   * @brief Default destructor.
   */
  ~Ellipse();

 private:
  // The center of the ellipse.
  const Eigen::Vector3d _center;
  // The first radius.
  Eigen::Vector3d _rad1;
  // The second radius.
  Eigen::Vector3d _rad2;
  // The plane in which the ellipse is.
  std::shared_ptr<Plane> _plane;
  // Get a rotation matrix around n by the angle a.
  Eigen::Matrix3d getRotationMatrixinR3(Eigen::Vector3d n, double a);
  // Integrate under the function of the ellipse (transformed to the origin, starting from -a)
  double integrate(double x, double a, double b);
  // Calculate the center of gravity for a fracture of the ellipse by integration.
  Eigen::Vector3d centerOfGravityEllipseFragment(double end, double a, double b, int ySign);
  // Integrate the x-coordinate for the center of gravity (ellipse is transformed to the origin).
  double integrateXY(double x, double a, double b);
  // Integrate the y-coordinate for the center of gravity (ellipse is transformed to the origin).
  double integrateYY(double x, double a, double b);
};

} /* namespace Serenity */

#endif /* GEOMETRY_ELLIPSE_H_ */
