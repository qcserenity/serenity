/**
 * @file   LBFGS.cpp
 *
 * @date   July 12, 2017
 * @author Kevin Klahr
 *  * @copyright \n
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
#include "math/optimizer/LBFGS.h"
/* Include Serenity Internal Headers */

/* Include Std and External Headers */

namespace Serenity {

LBFGS::LBFGS(Eigen::VectorXd& parameters) : Optimizer(parameters) {
}

void LBFGS::optimize(
    std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
    std::shared_ptr<unsigned int> nRejected) {
  (void)nRejected;
  /* maxm might need some tweaking */
  constexpr unsigned int maxm = 50;

  /* steplength for initial steepest descent step */
  double stepLength = 1.0;
  /* number of parameters treated */
  unsigned int nParams = this->_parameters.size();
  double value;
  unsigned int m = 0;
  bool converged = false;

  Eigen::MatrixXd dg(nParams, maxm);
  dg.setZero();
  Eigen::MatrixXd dx(nParams, maxm);
  dx.setZero();

  /* Setting all old gradients to zero. */
  Eigen::VectorXd gradients(nParams);
  gradients.setZero();
  /* Setting all old gradients to zero. */
  Eigen::VectorXd gradientsOld(nParams);
  gradientsOld.setZero();
  /* Set "assumed previous parameters" as parametersOld. */
  Eigen::VectorXd parametersOld(nParams);
  parametersOld = this->_parameters;

  converged = updateFunction(this->_parameters, value, gradients, nullptr, true);

  Eigen::VectorXd stepVector = 0.1 * gradients;
  while (true) {
    if (converged)
      break;
    /* Assign the (now) old parameter values to the private variables. */
    parametersOld = this->_parameters;
    /* Assign the old gradient values to the private variables. */
    gradientsOld = gradients;
    this->_parameters -= stepLength * stepVector;
    converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
    stepVector = gradients;
    if (m < maxm) {
      dg.col(m) = gradients - gradientsOld;
      dx.col(m) = this->_parameters - parametersOld;
      ++m;
    }
    else {
      dg.leftCols(maxm - 1) = dg.rightCols(maxm - 1);
      dx.leftCols(maxm - 1) = dx.rightCols(maxm - 1);
      dg.col(maxm - 1) = gradients - gradientsOld;
      dx.col(maxm - 1) = this->_parameters - parametersOld;
    }
    // L-BFGS update
    Eigen::VectorXd alpha(m);
    for (int i = m - 1; i > -1; --i) {
      double dxDotdg = dx.col(i).dot(dg.col(i));
      if (dxDotdg < 1.0e-3) {
        alpha[i] = dx.col(i).dot(stepVector) / 1.0e-3;
      }
      else {
        alpha[i] = dx.col(i).dot(stepVector) / dxDotdg;
      }
      stepVector -= alpha[i] * dg.col(i);
    }
    stepVector *= dx.col(m - 1).dot(dg.col(m - 1)) / dg.col(m - 1).dot(dg.col(m - 1));
    for (unsigned int i = 0; i < m; ++i) {
      double dxDotdg = dx.col(i).dot(dg.col(i));
      double beta = dg.col(i).dot(stepVector);
      if (dxDotdg < 1.0e-3) {
        beta /= 1.0e-3;
      }
      else {
        beta /= dx.col(i).dot(dg.col(i));
      }
      stepVector += dx.col(i) * (alpha[i] - beta);
    }
  }
}

} /* namespace Serenity */
