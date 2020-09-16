/**
 * @file   StaticDamping_test.cpp
 *
 * @date   Oct 14, 2016
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
#include "scf/damper/StaticDamping.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

namespace Serenity {
using namespace std;

/**
 * @test
 * @brief Tests StaticDamping.h: Restricted
 */
TEST(StaticDampingTest, Restricted_StaticDamping) {
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::MatrixXd> test(3, 3);

  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(2, 0) = 3.0;
  test(0, 1) = 4.0;
  test(1, 1) = 5.0;
  test(2, 1) = 6.0;
  test(0, 2) = 7.0;
  test(1, 2) = 8.0;
  test(2, 2) = 9.0;
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::MatrixXd> test2(3, 3);
  test2(0, 0) = 9.0;
  test2(1, 0) = 8.0;
  test2(2, 0) = 7.0;
  test2(0, 1) = 6.0;
  test2(1, 1) = 5.0;
  test2(2, 1) = 4.0;
  test2(0, 2) = 3.0;
  test2(1, 2) = 2.0;
  test2(2, 2) = 1.0;

  StaticDamping<Options::SCF_MODES::RESTRICTED> damper(0.5);

  damper.damp(test);
  damper.damp(test2);

  EXPECT_EQ(test2(0, 0), 5.0);
  EXPECT_EQ(test2(1, 0), 5.0);
  EXPECT_EQ(test2(2, 0), 5.0);
  EXPECT_EQ(test2(0, 1), 5.0);
  EXPECT_EQ(test2(1, 1), 5.0);
  EXPECT_EQ(test2(2, 1), 5.0);
  EXPECT_EQ(test2(0, 2), 5.0);
  EXPECT_EQ(test2(1, 2), 5.0);
  EXPECT_EQ(test2(2, 2), 5.0);
}

/**
 * @test
 * @brief Tests StaticDamping.h: Unrestricted
 */
TEST(StaticDampingTest, Unrestricted_StaticDamping) {
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> test(3, 3);
  test.alpha(0, 0) = 1.0;
  test.alpha(1, 0) = 2.0;
  test.alpha(2, 0) = 3.0;
  test.alpha(0, 1) = 4.0;
  test.alpha(1, 1) = 5.0;
  test.alpha(2, 1) = 6.0;
  test.alpha(0, 2) = 7.0;
  test.alpha(1, 2) = 8.0;
  test.alpha(2, 2) = 9.0;

  test.beta(0, 0) = 1.0;
  test.beta(1, 0) = 2.0;
  test.beta(2, 0) = 3.0;
  test.beta(0, 1) = 4.0;
  test.beta(1, 1) = 5.0;
  test.beta(2, 1) = 6.0;
  test.beta(0, 2) = 7.0;
  test.beta(1, 2) = 8.0;
  test.beta(2, 2) = 9.0;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> test2(3, 3);
  test2.alpha(0, 0) = 9.0;
  test2.alpha(1, 0) = 8.0;
  test2.alpha(2, 0) = 7.0;
  test2.alpha(0, 1) = 6.0;
  test2.alpha(1, 1) = 5.0;
  test2.alpha(2, 1) = 4.0;
  test2.alpha(0, 2) = 3.0;
  test2.alpha(1, 2) = 2.0;
  test2.alpha(2, 2) = 1.0;

  test2.beta(0, 0) = 9.0;
  test2.beta(1, 0) = 8.0;
  test2.beta(2, 0) = 7.0;
  test2.beta(0, 1) = 6.0;
  test2.beta(1, 1) = 5.0;
  test2.beta(2, 1) = 4.0;
  test2.beta(0, 2) = 3.0;
  test2.beta(1, 2) = 2.0;
  test2.beta(2, 2) = 1.0;

  StaticDamping<Options::SCF_MODES::UNRESTRICTED> damper(0.5);

  damper.damp(test);
  damper.damp(test2);

  EXPECT_EQ(test2.alpha(0, 0), 5.0);
  EXPECT_EQ(test2.alpha(1, 0), 5.0);
  EXPECT_EQ(test2.alpha(2, 0), 5.0);
  EXPECT_EQ(test2.alpha(0, 1), 5.0);
  EXPECT_EQ(test2.alpha(1, 1), 5.0);
  EXPECT_EQ(test2.alpha(2, 1), 5.0);
  EXPECT_EQ(test2.alpha(0, 2), 5.0);
  EXPECT_EQ(test2.alpha(1, 2), 5.0);
  EXPECT_EQ(test2.alpha(2, 2), 5.0);

  EXPECT_EQ(test2.beta(0, 0), 5.0);
  EXPECT_EQ(test2.beta(1, 0), 5.0);
  EXPECT_EQ(test2.beta(2, 0), 5.0);
  EXPECT_EQ(test2.beta(0, 1), 5.0);
  EXPECT_EQ(test2.beta(1, 1), 5.0);
  EXPECT_EQ(test2.beta(2, 1), 5.0);
  EXPECT_EQ(test2.beta(0, 2), 5.0);
  EXPECT_EQ(test2.beta(1, 2), 5.0);
  EXPECT_EQ(test2.beta(2, 2), 5.0);
}

} /* namespace Serenity */
