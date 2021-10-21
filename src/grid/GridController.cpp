/**
 * @file   GridController.cpp
 *
 * @date   Apr 29, 2014
 * @author Thomas Dresselhaus
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
#include "grid/GridController.h"

/* Further includes */

namespace Serenity {

GridController::GridController(std::unique_ptr<Grid> grid) : _grid(std::move(grid)) {
  assert(_grid);
}

const Eigen::Matrix3Xd& GridController::getGridPoints() {
  if (!_grid) {
    produceGrid();
  }
  return _grid->_gridPoints;
}

const Eigen::VectorXd& GridController::getWeights() {
  if (!_grid) {
    produceGrid();
  }
  return _grid->_weights;
}

unsigned int GridController::getNGridPoints() {
  if (!_grid) {
    produceGrid();
  }
  return _grid->_weights.size();
}

} /* namespace Serenity */
