/**
 * @file BFGS.cpp
 *
 * @date Oct 12, 2015
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
#include "math/optimizer/BFGS.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h" //Filtered output.
/* Include Std and External Headers */
#include <algorithm>
#include <iostream>

namespace Serenity {

BFGS::BFGS(Eigen::VectorXd& parameters, double initialStepLength, bool lineSearch, Eigen::MatrixXd initialInverseHess)
  : Optimizer(parameters), _initialStepLength(initialStepLength), _lineSearch(lineSearch), _initialInverseHess(initialInverseHess) {
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter" // value is ignored -> not needed
void BFGS::optimize(
    std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
    std::shared_ptr<unsigned int> nRejected) {
  (void)nRejected;
  /*
   * Get first set of parameters, gradients
   */
  double value = std::numeric_limits<double>::infinity();
  bool converged = false;
  /*
   * Setup variables
   */
  unsigned int nParams = this->_parameters.size();
  /*
   * Set hessInverse to identity -> First step is steepest descent
   */
  Eigen::VectorXd gradients(nParams);
  gradients.setZero();
  Eigen::MatrixXd hessInverse(nParams, nParams);
  if (_initialInverseHess.rows() == 0) {
    hessInverse.setIdentity();
  }
  else {
    hessInverse = _initialInverseHess;
  }
  Eigen::VectorXd gradientsOld(nParams);
  gradientsOld.setZero();
  Eigen::VectorXd parametersOld(nParams);
  parametersOld.setZero();
  Eigen::MatrixXd identity(nParams, nParams);
  identity.setIdentity();
  double oldValue = 0.0;
  double stepLength = _initialStepLength;
  // Estimated values for c1 and c2 in the Wolfe Conditions.
  // I have found multiple accounts claiming c1 should be very small (around 0.0001)
  // while some wrote it should be close to 0.5 (With the original publication putting it at 0.5).
  // For c2, no conclusive suggestion was found, 0.9 is a guess
  double c1 = 0.0001;
  double c2 = 0.9;

  converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
  /*
   * Check gradients whether their lengths match.
   */
  assert(nParams == gradients.size());
  /*
   * Update gradients, values
   */
  while (true) {
    if (converged) {
      break;
    }
    /*
     * Set old variables
     */
    parametersOld = this->_parameters;
    oldValue = value;
    /*
     * Build stepvector and update parameters
     */
    Eigen::VectorXd stepVector = -1.0 * hessInverse * gradients;
    // Some optimizations run into problems with stepsizes updating over multiple cycles
    if (!_lineSearch) {
      stepLength = 1.0;
    }
    if (!_lineSearch && stepVector.norm() > 0.5) {
      stepLength = 0.5 / stepVector.norm();
    }
    this->_parameters += stepLength * stepVector;
    gradientsOld = gradients;
    converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
    if (_lineSearch) {
      /*
       * Armijo Condition, keep steplength sufficiently small
       */
      if (value <= (oldValue + c1 * stepLength * gradientsOld.dot(stepVector))) {
        /*
         * Curvature Condition, prevent steplength from becoming too small
         */
        if (gradients.dot(stepVector) >= (c2 * gradientsOld.dot(stepVector))) {
          /*
           * Build new hessInverse by approximation through the Sherman-Morrison formula
           */
          Eigen::VectorXd parametersDiff = this->_parameters - parametersOld;
          Eigen::VectorXd gradientsDiff = gradients - gradientsOld;
          double paramDotGradInverse = 1 / (parametersDiff.dot(gradientsDiff));
          if ((parametersDiff.dot(gradientsDiff)) < 1e-6) {
            paramDotGradInverse = 100000;
          }
          Eigen::MatrixXd A1 = identity - (parametersDiff * gradientsDiff.transpose()) * paramDotGradInverse;
          Eigen::MatrixXd A2 = identity - (gradientsDiff * parametersDiff.transpose()) * paramDotGradInverse;
          hessInverse = A1 * (hessInverse * A2) + paramDotGradInverse * parametersDiff * parametersDiff.transpose();
        }
        else {
          // Backtrack if energy increased, otherwise just keep the step to save time despite requiring a change in the
          // steplength
          if (value > oldValue) {
            OutputControl::vOut << "Energy increased... Backtracking..." << std::endl;
            this->_parameters -= stepLength * stepVector;
            gradients = gradientsOld;
          }
          else {
            /*
             * Build new hessInverse by approximation through the Sherman-Morrison formula
             */
            Eigen::VectorXd parametersDiff = this->_parameters - parametersOld;
            Eigen::VectorXd gradientsDiff = gradients - gradientsOld;
            double paramDotGradInverse = 1 / (parametersDiff.dot(gradientsDiff));
            if ((parametersDiff.dot(gradientsDiff)) < 1e-6) {
              paramDotGradInverse = 100000;
            }
            Eigen::MatrixXd A1 = identity - (parametersDiff * gradientsDiff.transpose()) * paramDotGradInverse;
            Eigen::MatrixXd A2 = identity - (gradientsDiff * parametersDiff.transpose()) * paramDotGradInverse;
            hessInverse = A1 * (hessInverse * A2) + paramDotGradInverse * parametersDiff * parametersDiff.transpose();
          }
          OutputControl::vOut << "Increasing step length" << std::endl;
          stepLength *= 1.5;
        }
      }
      else {
        // Backtrack if energy increased, otherwise just keep the step to save time despite requiring a change in the
        // steplength
        if (value > oldValue) {
          OutputControl::vOut << "Energy increased... Backtracking..." << std::endl;
          this->_parameters -= stepLength * stepVector;
          gradients = gradientsOld;
        }
        else {
          /*
           * Build new hessInverse by approximation through the Sherman-Morrison formula
           */
          Eigen::VectorXd parametersDiff = this->_parameters - parametersOld;
          Eigen::VectorXd gradientsDiff = gradients - gradientsOld;
          double paramDotGradInverse = 1 / (parametersDiff.dot(gradientsDiff));
          if ((parametersDiff.dot(gradientsDiff)) < 1e-6) {
            paramDotGradInverse = 100000;
          }
          Eigen::MatrixXd A1 = identity - (parametersDiff * gradientsDiff.transpose()) * paramDotGradInverse;
          Eigen::MatrixXd A2 = identity - (gradientsDiff * parametersDiff.transpose()) * paramDotGradInverse;
          hessInverse = A1 * (hessInverse * A2) + paramDotGradInverse * parametersDiff * parametersDiff.transpose();
        }
        OutputControl::vOut << "Decreasing step length" << std::endl;
        stepLength *= 0.5;
      }
    }
    else {
      /*
       * Build new hessInverse by approximation through the Sherman-Morrison formula
       */
      Eigen::VectorXd parametersDiff = this->_parameters - parametersOld;
      Eigen::VectorXd gradientsDiff = gradients - gradientsOld;
      double paramDotGradInverse = 1 / (parametersDiff.dot(gradientsDiff));
      if ((parametersDiff.dot(gradientsDiff)) < 1e-6) {
        paramDotGradInverse = 100000;
      }
      Eigen::MatrixXd A1 = identity - (parametersDiff * gradientsDiff.transpose()) * paramDotGradInverse;
      Eigen::MatrixXd A2 = identity - (gradientsDiff * parametersDiff.transpose()) * paramDotGradInverse;
      hessInverse = A1 * (hessInverse * A2) + paramDotGradInverse * parametersDiff * parametersDiff.transpose();
    }
  }
}

#pragma GCC diagnostic pop // "-Wunused-parameter"

} /* namespace Serenity */