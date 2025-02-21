/**
 * @file   Optimizer.h
 *
 * @date   Jul 21, 2015
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
#ifndef MATHS_OPTIMIZER_OPTIMIZER_H_
#define MATHS_OPTIMIZER_OPTIMIZER_H_
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <cassert>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
/**
 * @brief Possible modes for an Optimizer, meaning minimization or maximization
 */
enum class OptimizationMode { MINIMIZE, MAXIMIZE };

/**
 * @class Optimizer Optimizer.h
 * @brief Base class for all optimization algorithms
 */
class Optimizer {
 public:
  /**
   * @brief Constructor
   * @param mode The mode of the optimization e.g. maximization, default: minimization
   */
  Optimizer(Eigen::VectorXd& parameters, OptimizationMode mode = OptimizationMode::MINIMIZE)
    : _parameters(parameters), _mode(mode), _allowModeSwitches(true) {
  }
  /**
   * @brief Default destructor.
   */
  virtual ~Optimizer() = default;
  /**
   * @brief Performs a full optimization of a set of parameters along a given
   *        updateFunction.
   * @param updateFunction Lambda function that updates energies and gradients at a given position.
   *                       The lambda function should use a set of parameters (1. argument) to evaluate
   *                       the value and gradient (2. and 3. argument) of the function.
   * @param nRejected Counter to keep track of the number of rejected steps and obsolete cycles.
   * A short example:
   * @code
   * auto updateFunction=[](const Eigen::VectorXd& parameters,
   *                        double& value,
   *                        Eigen::VectorXd& gradients)->void{
   *   Update value and gradients here!
   * };
   * Eigen::vectorXd parameters;
   * double value;
   * Eigen::VectorXd gradients;
   *
   * Optimizer opt(parameters, value, gradients);
   * opt.optimize(updateFunction);
   * @endcode
   */
  virtual void optimize(std::function<bool(const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                           std::shared_ptr<Eigen::MatrixXd> hessian, bool print)>
                            updateFunction,
                        std::shared_ptr<unsigned int> nRejected = nullptr) = 0;
  /**
   * @brief Sets the mode of the optimization.
   * @param mode The new mode for the optimizer.
   *
   * Make sure that either the optimization mode is actually interchangeable at all times,
   * or disable it with the _allowModeSwitches boolean. This can be done inside the constructor
   * or (better) inside the optimize routine.
   */
  void setOptimizationMode(OptimizationMode mode) {
    assert(_allowModeSwitches);
    _mode = mode;
  };

  virtual void reinit(){};

 protected:
  Eigen::VectorXd& _parameters;
  OptimizationMode _mode;
  bool _allowModeSwitches;
};

} /* namespace Serenity */

#endif /* MATHS_OPTIMIZER_OPTIMIZER_H_ */
