/**
 * @file   GridController__TEST_SUPPLY.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   23. August 2015, 18:26
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
/* Include Class Header*/
#include "testsupply/GridController__TEST_SUPPLY.h"
/* Include Serenity Internal Headers */
#include "grid/Grid.h"

namespace Serenity {

std::map<TEST_GRID_CONTROLLERS, std::shared_ptr<GridController>> GridController__TEST_SUPPLY::_testGridControllers;

void GridController__TEST_SUPPLY::prepare(TEST_GRID_CONTROLLERS kind) {
  switch (kind) {
    case TEST_GRID_CONTROLLERS::TINY: {
      std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, 4));
      (*points) << 0.0, 1.0, -2.0, 3.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0, 1.0, -1.0;
      std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(4));
      (*weights) << 1.0, 1.0, 1.0, 1.0;
      _testGridControllers[kind] =
          std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
    } break;
    case TEST_GRID_CONTROLLERS::VERY_SMALL: {
      std::unique_ptr<Eigen::Matrix3Xd> points(new Eigen::Matrix3Xd(3, 5));
      (*points) << -1.0, 0.0, 0.0, 0.8, 0.1, 0.0, 0.5, -0.4, 0.1, -0.1, 0.0, 0.4, -0.5, -0.1, 0.1;
      std::unique_ptr<Eigen::VectorXd> weights(new Eigen::VectorXd(5));
      (*weights) << 0.1, 0.2, 0.15, 0.15, 0.4;
      _testGridControllers[kind] =
          std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(points), std::move(weights))));
    } break;
  }
}

} /* namespace Serenity */
