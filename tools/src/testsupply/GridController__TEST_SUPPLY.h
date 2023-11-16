/**
 * @file   GridController__TEST_SUPPLY.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   23. August 2015, 17:32
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
#ifndef GRIDCONTROLLER__TEST_SUPPLY_H
#define GRIDCONTROLLER__TEST_SUPPLY_H
/* Include Serenity Internal Headers */
#include "grid/GridController.h"
/* Include Std and External Headers */
#include <map>
#include <memory>

namespace Serenity {
/**
 * These kinds of test grids are available:\n
 *
 * TINY:\n
 *
 * Not a proper integration grid. Just 4 points all with weight 1.
 *
 * VERY_SMALL:\n
 *
 * Also not a proper integration grid. But weights sum up to 1.0. 5 points.
 */
enum class TEST_GRID_CONTROLLERS { TINY, VERY_SMALL };
/**
 * @class GridController__TEST_SUPPLY GridController__TEST_SUPPLY.h
 * @brief Provides objects of type GridController without further dependencies.
 *
 * I.e. The grids are ready for usage.
 */
class GridController__TEST_SUPPLY {
  GridController__TEST_SUPPLY() = delete;

 public:
  static std::shared_ptr<GridController> getGridController(TEST_GRID_CONTROLLERS kind) {
    if (!_testGridControllers[kind])
      prepare(kind);
    return _testGridControllers[kind];
  }

 private:
  static void prepare(TEST_GRID_CONTROLLERS kind);
  static std::map<TEST_GRID_CONTROLLERS, std::shared_ptr<GridController>> _testGridControllers;
};

} /* namespace Serenity */

#endif /* GRIDCONTROLLER__TEST_SUPPLY_H */
