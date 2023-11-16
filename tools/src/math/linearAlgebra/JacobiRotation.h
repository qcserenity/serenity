/**
 * @file JacobiRotation.h
 *
 * @date Mar 10, 2019
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

#ifndef MATH_LINEARALGEBRA_JACOBIROTATION_H_
#define MATH_LINEARALGEBRA_JACOBIROTATION_H_

/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {
/**
 * @class JacobiRotation JacobiRotation.h
 *
 * @brief A class to perform single 2x2 Jacobi rotation.
 */
class JacobiRotation {
 public:
  /**
   * @brief Constructor
   */
  JacobiRotation() = default;
  /**
   * @brief Destructor
   */
  virtual ~JacobiRotation() = default;

  /**
   * @brief Performs a single 2x2 Jacobi rotation.
   * @param vec1 First vector to be rotated.
   * @param vec2 Second vector to be rotated.
   * @param angle The rotation angle.
   */
  static void rotate(Eigen::Ref<Eigen::ArrayXd> vec1, Eigen::Ref<Eigen::ArrayXd> vec2, double angle);
};

} /* namespace Serenity */

#endif /* MATH_LINEARALGEBRA_JACOBIROTATION_H_ */
