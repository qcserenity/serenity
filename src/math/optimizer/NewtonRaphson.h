/**
 * @file   NewtonRaphson.h
 *
 * @date   Jul 22, 2015
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef MATHS_OPTIMIZER_NEWTONRAPHSON_H_
#define MATHS_OPTIMIZER_NEWTONRAPHSON_H_
/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"


namespace Serenity {
/* Forward declarations */
class Lapack;
template<class T>class Matrix;
/**
 * @class NewtonRaphson NewtonRaphson
 * @brief Derived Optimizer class for the Newton-Raphson algorithm.
 *
 * This class can be used to run Newton-Raphson optimizations.
 * because this type of algorithm is uses a Hessian for each update step it has one
 * major set of additional functions.
 * The main idea here is to keep the Newthon-Raphson Optimizer as a derived class of
 * Optimizer while not requiring every other Optimizer to be handed a Hessian in
 * every step.
 *
 */
class NewtonRaphson : public Optimizer {
public:
  NewtonRaphson(Eigen::VectorXd& parameters, double threshold=0.0);
  virtual ~NewtonRaphson() = default;

  /**
   * @brief See documentation of base class Optimizer.
   */
  void optimize(std::function<bool(
      const Eigen::VectorXd&,
      double&,
      Eigen::VectorXd&,
      std::shared_ptr<Eigen::MatrixXd> hessian,
      bool print)>updateFunction) override final;

  /**
   * @brief Set the threshold for the minimal eigenvalue before a pseudoinversion is performed
   * @param threshold The threshold for the eigenvalues befor a pseudoinversion is performed
   */

  void setThreshold(double threshold){
    _threshold=threshold;
  }

private:
  double _threshold;
};

} /* namespace Serenity */

#endif /* MATHS_OPTIMIZER_NEWTONRAPHSON_H_ */
