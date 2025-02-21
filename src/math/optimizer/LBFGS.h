/**
 * @file   LBFGS.h
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
#ifndef MATH_OPTIMIZER_LBFGS_H_
#define MATH_OPTIMIZER_LBFGS_H_

/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"
/* Include Std and External Headers */

namespace Serenity {

/**
 * @class LBFGS LBFGS.h
 * @brief Implementation of the Low-memory Broyden–Fletcher–Goldfarb–Shanno (BFGS) optimization algorithm.
 */
class LBFGS : public Optimizer {
 public:
  LBFGS(Eigen::VectorXd& parameters);
  virtual ~LBFGS() = default;

  /**
   * @brief See documentation of base class Optimizer.
   *        Note that no Hessian matrix is required and a nullptr is passed to the optimize function.
   */
  void
  optimize(std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
           std::shared_ptr<unsigned int> nRejected = nullptr) override final;
};
} /* namespace Serenity */

#endif /* MATH_OPTIMIZER_LBFGS_H_ */
