/**
 * @file JacobiRotation.cpp
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

/* Include Class Header*/
#include "math/linearAlgebra/JacobiRotation.h"

namespace Serenity {

void JacobiRotation::rotate(Eigen::Ref<Eigen::ArrayXd> vec1, Eigen::Ref<Eigen::ArrayXd> vec2, double angle) {
  const double cosAngle = cos(angle);
  const double sinAngle = sin(angle);
  Eigen::VectorXd bak = vec1;
  vec1.array() = cosAngle * bak.array() + sinAngle * vec2.array();
  vec2.array() = -sinAngle * bak.array() + cosAngle * vec2.array();
}

} /* namespace Serenity */
