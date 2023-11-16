/**
 * @file   HilbertRTreeSorting_test.cpp
 *
 * @date    Sep 29, 2017
 * @author  Jan Unsleber
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
#include "grid/HilbertRTreeSorting.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class HilbertRTreeSortingTest : public ::testing::Test {
 protected:
};

TEST_F(HilbertRTreeSortingTest, 2x2x2) {
  Eigen::VectorXd weights(8);
  weights << 0, 1, 2, 3, 4, 5, 6, 7;
  Eigen::Matrix3Xd coords(3, 8);
  coords << +0.5, -0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, +0.5, -0.5, +0.5, -0.5,
      -0.5, +0.5, -0.5, +0.5, -0.5, +0.5;

  HilbertRTreeSorting sort(coords, weights);
  sort.sort();

  Eigen::VectorXd wref(8);
  wref << 3, 7, 1, 4, 6, 2, 5, 0;
  Eigen::Matrix3Xd pref(3, 8);
  pref << -0.5, -0.5, -0.5, -0.5, +0.5, +0.5, +0.5, +0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, -0.5, +0.5, +0.5, +0.5,
      -0.5, -0.5, -0.5, -0.5, +0.5, +0.5;

  for (unsigned int i = 0; i < 8; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      EXPECT_EQ(pref(j, i), coords(j, i));
    }
    EXPECT_EQ(wref[i], weights[i]);
  }
}

TEST_F(HilbertRTreeSortingTest, 2x2x2_2cores) {
  Eigen::VectorXd weights(8);
  weights << 0, 1, 2, 3, 4, 5, 6, 7;
  Eigen::Matrix3Xd coords(3, 8);
  coords << +0.5, -0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, +0.5, -0.5, +0.5, -0.5,
      -0.5, +0.5, -0.5, +0.5, -0.5, +0.5;

  int threads = omp_get_max_threads();
  omp_set_num_threads(2);
  HilbertRTreeSorting sort(coords, weights);
  sort.sort();
  omp_set_num_threads(threads);

  Eigen::VectorXd wref(8);
  wref << 3, 7, 1, 4, 6, 2, 5, 0;
  Eigen::Matrix3Xd pref(3, 8);
  pref << -0.5, -0.5, -0.5, -0.5, +0.5, +0.5, +0.5, +0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5, -0.5, +0.5, +0.5, +0.5,
      -0.5, -0.5, -0.5, -0.5, +0.5, +0.5;

  for (unsigned int i = 0; i < 8; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      EXPECT_EQ(pref(j, i), coords(j, i));
    }
    EXPECT_EQ(wref[i], weights[i]);
  }
}

} /* namespace Serenity */
