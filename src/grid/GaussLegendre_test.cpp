/**
 * @file   GaussLegendre_test.cpp
 *
 * @date   10.04.2020
 * @author J. Toelle
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

/* Include Serenity Internal Headers */
#include "grid/GaussLegendre.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

namespace Serenity {
TEST(GaussLegendreTest, Integration) {
  auto integrator = GaussLegendre(128);
  auto weights = integrator.getWeights();
  auto nodes = integrator.getGridPoints();

  // Integral transformation for integration from [-inf,inf]
  Eigen::VectorXd nodes_Modified = (1.0 + nodes.array()) / (1.0 - nodes.array());
  weights = (2.0 * weights.array()) / ((1.0 - nodes.array()).array().square());

  double integrant = 0.0;
  // Integration from [0,inf]
  for (unsigned int i = 0; i < weights.size(); i++) {
    integrant += weights(i) * (1.0 / (pow(nodes_Modified(i), 2) + 0.1));
  }
  integrant = 2.0 * integrant;
  std::cout << " The integrant evaluates to -- " << integrant << std::endl;
  // Value from wolfrahm alpha (10.04.2020)
  EXPECT_NEAR(integrant, 9.9345882657961012344335506703, 1e-14);
}
} /* namespace Serenity */