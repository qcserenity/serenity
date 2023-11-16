/**
 * @file   BFGS_test.cpp
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
/* Include Serenity Internal Headers */
#include "math/optimizer/BFGS.h"
#include "math/optimizer/LBFGS.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

namespace Serenity {

TEST(BFGSTest, BFGSConvergence) {
  auto f = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return vec[0] * vec[0] + 2 * vec[1] * vec[1] - 0.3 * cos(3 * pi * vec[0]) - 0.4 * cos(4 * pi * vec[1]) + 0.7;
  };
  auto dfdx = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return 2 * vec[0] + 0.3 * sin(3 * pi * vec[0]);
  };
  auto dfdy = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return 4 * vec[1] * vec[1] + 0.4 * sin(4 * pi * vec[1]);
  };

  Eigen::VectorXd vec(2);
  vec << 0.5, 1.0;

  BFGS optimizer(vec);

  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
    (void)hessian;
    (void)printInfo;
    gradients[0] = dfdx(parameters);
    gradients[1] = dfdy(parameters);
    value = f(parameters);

    bool converged = false;
    if (gradients.norm() < 1e-6) {
      converged = true;
    }

    return converged;
  };

  optimizer.optimize(updateFunction);

  EXPECT_NEAR(vec[0], 0, 1e-6);
  EXPECT_NEAR(vec[1], 0, 1e-6);
}

TEST(LBFGSTest, LBFGSConvergence) {
  auto f = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return vec[0] * vec[0] + 2 * vec[1] * vec[1] - 0.3 * cos(3 * pi * vec[0]) - 0.4 * cos(4 * pi * vec[1]) + 0.7;
  };
  auto dfdx = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return 2 * vec[0] + 0.3 * sin(3 * pi * vec[0]);
  };
  auto dfdy = [](const Eigen::VectorXd& vec) {
    double pi = 3.1415926535;
    return 4 * vec[1] * vec[1] + 0.4 * sin(4 * pi * vec[1]);
  };

  Eigen::VectorXd vec(2);
  vec << 0.5, 1.0;

  LBFGS optimizer(vec);

  unsigned int cycle = 0;
  auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                  std::shared_ptr<Eigen::MatrixXd> hessian, bool printInfo) {
    (void)hessian;
    (void)printInfo;
    gradients[0] = dfdx(parameters);
    gradients[1] = dfdy(parameters);
    value = f(parameters);

    bool converged = false;
    if (gradients.norm() < 1e-6) {
      converged = true;
    }

    cycle++;

    return converged;
  };

  optimizer.optimize(updateFunction);

  EXPECT_NEAR(vec[0], 0, 1e-6);
  EXPECT_NEAR(vec[1], 0, 1e-6);
}

} /* namespace Serenity */
