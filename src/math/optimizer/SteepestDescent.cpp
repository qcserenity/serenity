/**
 * @file   SteepestDescent.cpp
 *
 * @date   Oct 12, 2015
 * @author Linus Scholz
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
#include "math/optimizer/SteepestDescent.h"

namespace Serenity {

SteepestDescent::SteepestDescent(Eigen::VectorXd& parameters, double stepWidth)
  : Optimizer(parameters), _stepWidth(stepWidth) {
}
void SteepestDescent::optimize(
    std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
    std::shared_ptr<unsigned int> nRejected) {
  (void)nRejected;
  double value = std::numeric_limits<double>::infinity();
  Eigen::VectorXd gradients(this->_parameters.size());
  gradients.setZero();
  bool converged = false;

  do {
    this->_parameters -= _stepWidth * gradients;
    converged = updateFunction(this->_parameters, value, gradients, nullptr, false);
  } while (!converged);
}
} /* namespace Serenity */
