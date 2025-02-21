/**
 * @file   SteepestDescent.h
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
#ifndef MATHS_OPTIMIZER_STEEPESTDESCENT_H_
#define MATHS_OPTIMIZER_STEEPESTDESCENT_H_
/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"

namespace Serenity {
/**
 * @class SteepestDescent SteepestDescent.h
 * @brief A simple steepest descent optimizer
 *
 * A steepest descent optimizer, the most simple optimizer possible
 * updates the parameters by a constant factor (step width) times the
 * gradients.
 *
 */
class SteepestDescent : public Optimizer {
 public:
  /**
   * @brief Constructor.
   * @param stepWidth The step width: the factor to multiply the gradients by.
   */
  SteepestDescent(Eigen::VectorXd& parameters, double stepWidth = 0.01);
  virtual ~SteepestDescent() = default;

  /**
   * @brief See documentation of base class Optimizer.
   */
  void
  optimize(std::function<bool(const Eigen::VectorXd&, double&, Eigen::VectorXd&, std::shared_ptr<Eigen::MatrixXd> hessian, bool print)> updateFunction,
           std::shared_ptr<unsigned int> nRejected = nullptr) override final;

 private:
  double _stepWidth;
};

} /* namespace Serenity */

#endif /* MATHS_OPTIMIZER_STEEPESTDESCENT_H_ */
