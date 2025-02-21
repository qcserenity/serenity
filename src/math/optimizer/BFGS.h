/**
 * @file BFGS.h
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
#ifndef MATHS_OPTIMIZER_BFGS_H_
#define MATHS_OPTIMIZER_BFGS_H_
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
#include "math/optimizer/Optimizer.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/**
 * @class BFGS BFGS.h
 * @brief Implementation of the Broyden–Fletcher–Goldfarb–Shanno (BFGS) optimization algorithm.
 */
class BFGS : public Optimizer {
 public:
  /**
   * @brief Constructor.
   * @param parameters The parameters to be optimized.
   * @param initialStepLength The initial step length.
   * @param lineSearch Flag to turn on a line search algorithm.
   * @param initialInverseHess Starting Hessian.
   */
  BFGS(Eigen::VectorXd& parameters, double initialStepLength = 1.0, bool lineSearch = false,
       Eigen::MatrixXd initialInverseHess = Eigen::MatrixXd(0, 0));
  /**
   * @brief Default destructor.
   */
  virtual ~BFGS() = default;

  /**
   * @brief See documentation of base class Optimizer.
   *
   * The update function parameter 'hessian' is ignored.
   */
  void optimize(std::function<bool(const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                   std::shared_ptr<Eigen::MatrixXd> hessian, bool print)>
                    updateFunction,
                std::shared_ptr<unsigned int> nRejected = nullptr) override final;

 private:
  // The initial step length.
  double _initialStepLength;
  // Flag for using a line search algorithm.
  bool _lineSearch;
  // The starting Hessian.
  Eigen::MatrixXd _initialInverseHess;
};

} /* namespace Serenity */

#endif /* MATHS_OPTIMIZER_BFGS_H_ */
