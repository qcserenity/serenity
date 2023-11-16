/**
 * @file NewtonRaphson_test.cpp
 *
 * @date Jul 22, 2015
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
/* Include Serenity Internal Headers */
#include "math/optimizer/NewtonRaphson.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests optimizer NewtonRaphson.h: for 4(x-1)**2 and x=5
 */
TEST(NewtonRaphsonTest, QuadraticConvergenceWrapper) {
  auto f = [](const Eigen::VectorXd& vec) {
    return (vec[0] * vec[0] * vec[0] * vec[0] - 16.0 * vec[0] * vec[0] + 5.0 * vec[0]) / 2.0 +
           (vec[1] * vec[1] * vec[1] * vec[1] - 16.0 * vec[1] * vec[1] + 5.0 * vec[1]) / 2.0;
  };
  auto dfdx = [](const Eigen::VectorXd& vec) {
    return (4.0 * vec[0] * vec[0] * vec[0] - 16.0 * 2.0 * vec[0] + 5.0) / 2.0;
  };
  auto dfdy = [](const Eigen::VectorXd& vec) {
    return (4.0 * vec[1] * vec[1] * vec[1] - 16.0 * 2.0 * vec[1] + 5.0) / 2.0;
  };
  auto d2fdx2 = [](const Eigen::VectorXd& vec) { return (4.0 * 3.0 * vec[0] * vec[0] - 16.0 * 2.0) / 2.0; };
  auto d2fdy2 = [](const Eigen::VectorXd& vec) { return (4.0 * 3.0 * vec[1] * vec[1] - 16.0 * 2.0) / 2.0; };

  Eigen::VectorXd vec(2);
  vec << 3.0, 3.2;

  NewtonRaphson optimizer(vec);

  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
    (void)printInfo;
    gradients[0] = dfdx(parameters);
    gradients[1] = dfdy(parameters);
    value = f(parameters);
    hessian->setZero();
    (*hessian)(0, 0) = d2fdx2(vec);
    (*hessian)(1, 1) = d2fdy2(vec);

    bool converged = false;
    if (gradients.norm() < 1e-4) {
      converged = true;
    }

    return converged;
  };

  optimizer.optimize(updateFunction);

  EXPECT_NEAR(vec[0], 2.7468027, 1e-6);
  EXPECT_NEAR(vec[1], 2.7468027, 1e-6);
}

} /* namespace Serenity */
