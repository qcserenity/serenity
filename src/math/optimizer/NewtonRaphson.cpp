/**
 * @file   NewtonRaphson.cpp
 *
 * @date   Jul 22, 2015
 * @author Jan Unsleber
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
#include "math/optimizer/NewtonRaphson.h"

namespace Serenity {

NewtonRaphson::NewtonRaphson(Eigen::VectorXd& parameters, double threshold)
  : Optimizer(parameters), _threshold(threshold) {
}

void NewtonRaphson::minimize2D(std::function<bool(const Eigen::VectorXd&, double&, double&, double& hessian, bool print)> updateFunction) {
  double value = std::numeric_limits<double>::infinity();
  bool converged = true;
  double gradient = 0.0;
  double hessian = 0.0;
  converged = updateFunction(this->_parameters, value, gradient, hessian, false);

  while (true) {
    if (converged)
      break;
    double steps = gradient / hessian;
    double sign = (hessian >= 0) ? 1.0 : -1.0;
    this->_parameters(0) -= sign * steps;
    converged = updateFunction(this->_parameters, value, gradient, hessian, false);
  }
}

void NewtonRaphson::optimize(
    std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
    std::shared_ptr<unsigned int> nRejected) {
  (void)nRejected;
  double value = std::numeric_limits<double>::infinity();
  bool converged = true;
  unsigned nParams = this->_parameters.size();
  Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
  auto hessian = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(nParams, nParams));
  converged = updateFunction(this->_parameters, value, gradients, hessian, false);

  while (true) {
    if (converged)
      break;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(*hessian, Eigen::ComputeFullV | Eigen::ComputeFullU);
    svd.setThreshold(_threshold);
    Eigen::VectorXd steps = svd.solve(gradients);
    this->_parameters -= steps;
    converged = updateFunction(this->_parameters, value, gradients, hessian, false);
  }
}

} /* namespace Serenity */
