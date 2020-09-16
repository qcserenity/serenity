/**
 * @file   QST_test.cpp
 *
 * @date   Jun 28, 2017
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
#include "math/saddlepoint/QST.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

namespace Serenity {

/**
 * @class QSTTest
 * @brief Sets everything up for the tests of QST.h .
 */
class QSTTest : public ::testing::Test {
 protected:
};

/**
 * @test
 * @brief Tests QST.h: LST guess then QST on an analytical function.
 */
TEST_F(QSTTest, LSTQST_Function) {
  Eigen::VectorXd mini1(2);
  mini1 << 0.6, -0.6;
  Eigen::VectorXd mini2(2);
  mini2 << 0.6, 1.1;

  double _value;
  Eigen::VectorXd _parameters;
  auto _getData = [&_value, &_parameters](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                          bool dontcare) -> void {
    (void)dontcare;
    assert(parameters.size() == 2);
    assert(gradients.size() == 2);
    Eigen::MatrixXd params(5, 8);
    params << 1.70, 1.70, 0.80, 0.80, -1.00, -0.25, -0.25, -0.50, 1.50, 1.50, 4.00, 4.00, 8.00, 2.00, 2.00, 4.00, -0.25,
        1.75, -0.10, 1.60, 0.55, -0.75, 2.25, 0.75, 4.00, 4.00, 4.00, 4.00, 2.00, 4.00, 4.00, 4.00, 0.50, 0.50, -0.95,
        -0.95, -0.75, -0.75, -0.75, 1.20;
    value = 0.0;
    gradients.setZero();
    for (int i = 0; i < 8; i++) {
      const double tmp =
          params(0, i) * exp(-params(1, i) * (parameters[0] - params(2, i)) * (parameters[0] - params(2, i)) -
                             params(3, i) * (parameters[1] - params(4, i)) * (parameters[1] - params(4, i)));
      value += tmp;
      gradients[0] += -2 * params(1, i) * (parameters[0] - params(2, i)) * tmp;
      gradients[1] += -2 * params(3, i) * (parameters[1] - params(4, i)) * tmp;
    }
    _value = value;
    _parameters = parameters;
  };

  QST qst(mini1, mini2, false);
  qst.optimize(_getData);
  EXPECT_NEAR(_parameters[0], 0.7279, 1e-04);
  EXPECT_NEAR(_parameters[1], 0.4733, 1e-04);
  auto t = qst.getTangent();
  EXPECT_NEAR(t[0], -0.084, 1e-3);
  EXPECT_NEAR(t[1], 0.996, 1e-3);
  EXPECT_NEAR(t.norm(), 1.0, 1e-7);
}
/**
 * @test
 * @brief Tests QST.h: User guess then QST on an analytical function.
 */
TEST_F(QSTTest, GuessQST_Function) {
  Eigen::VectorXd mini1(2);
  mini1 << 0.6, -0.6;
  Eigen::VectorXd mini2(2);
  mini2 << 0.6, 1.1;

  double _value;
  Eigen::VectorXd _parameters;
  auto _getData = [&_value, &_parameters](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                          bool dontcare) -> void {
    (void)dontcare;
    assert(parameters.size() == 2);
    assert(gradients.size() == 2);
    Eigen::MatrixXd params(5, 8);
    params << 1.70, 1.70, 0.80, 0.80, -1.00, -0.25, -0.25, -0.50, 1.50, 1.50, 4.00, 4.00, 8.00, 2.00, 2.00, 4.00, -0.25,
        1.75, -0.10, 1.60, 0.55, -0.75, 2.25, 0.75, 4.00, 4.00, 4.00, 4.00, 2.00, 4.00, 4.00, 4.00, 0.50, 0.50, -0.95,
        -0.95, -0.75, -0.75, -0.75, 1.20;
    value = 0.0;
    gradients.setZero();
    for (int i = 0; i < 8; i++) {
      const double tmp =
          params(0, i) * exp(-params(1, i) * (parameters[0] - params(2, i)) * (parameters[0] - params(2, i)) -
                             params(3, i) * (parameters[1] - params(4, i)) * (parameters[1] - params(4, i)));
      value += tmp;
      gradients[0] += -2 * params(1, i) * (parameters[0] - params(2, i)) * tmp;
      gradients[1] += -2 * params(3, i) * (parameters[1] - params(4, i)) * tmp;
    }
    _value = value;
    _parameters = parameters;
  };

  std::unique_ptr<Eigen::VectorXd> guess(new Eigen::VectorXd(2));
  (*guess)[0] = 0.0;
  (*guess)[1] = 0.5;

  QST qst(mini1, mini2, false, std::move(guess));
  qst.optimize(_getData);
  EXPECT_NEAR(_parameters[0], 0.7279, 1e-03);
  EXPECT_NEAR(_parameters[1], 0.4733, 1e-03);
}

/**
 * @test
 * @brief Tests QST.h: User guess then QST with the combined flag true, on an analytical function.
 */
TEST_F(QSTTest, GuessQST_Function_FewLoops) {
  Eigen::VectorXd mini1(2);
  mini1 << 0.6, -0.6;
  Eigen::VectorXd mini2(2);
  mini2 << 0.6, 1.1;

  double _value;
  Eigen::VectorXd _parameters;
  auto _getData = [&_value, &_parameters](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                          bool dontcare) -> void {
    (void)dontcare;
    assert(parameters.size() == 2);
    assert(gradients.size() == 2);
    Eigen::MatrixXd params(5, 8);
    params << 1.70, 1.70, 0.80, 0.80, -1.00, -0.25, -0.25, -0.50, 1.50, 1.50, 4.00, 4.00, 8.00, 2.00, 2.00, 4.00, -0.25,
        1.75, -0.10, 1.60, 0.55, -0.75, 2.25, 0.75, 4.00, 4.00, 4.00, 4.00, 2.00, 4.00, 4.00, 4.00, 0.50, 0.50, -0.95,
        -0.95, -0.75, -0.75, -0.75, 1.20;
    value = 0.0;
    gradients.setZero();
    for (int i = 0; i < 8; i++) {
      const double tmp =
          params(0, i) * exp(-params(1, i) * (parameters[0] - params(2, i)) * (parameters[0] - params(2, i)) -
                             params(3, i) * (parameters[1] - params(4, i)) * (parameters[1] - params(4, i)));
      value += tmp;
      gradients[0] += -2 * params(1, i) * (parameters[0] - params(2, i)) * tmp;
      gradients[1] += -2 * params(3, i) * (parameters[1] - params(4, i)) * tmp;
    }
    _value = value;
    _parameters = parameters;
  };

  std::unique_ptr<Eigen::VectorXd> guess(new Eigen::VectorXd(2));
  (*guess)[0] = 0.0;
  (*guess)[1] = 0.5;

  QST qst(mini1, mini2, true, std::move(guess));
  qst.optimize(_getData);
  EXPECT_NEAR(_parameters[0], 0.78618359296940343, 1e-04);
  EXPECT_NEAR(_parameters[1], 0.5301824320223022, 1e-04);
}

} /* namespace Serenity */
