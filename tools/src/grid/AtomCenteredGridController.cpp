/**
 * @file   AtomCenteredGridController.cpp
 *
 * @date   Mar 09, 2016
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
#include "grid/AtomCenteredGridController.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "grid/construction/GridFactory.h"

namespace Serenity {

AtomCenteredGridController::AtomCenteredGridController(const std::shared_ptr<const Geometry> geometry,
                                                       const std::shared_ptr<GridFactory> gridFactory,
                                                       const bool gridPointSorting)
  : _geometry(geometry), _gridFactory(gridFactory), _gridPointSorting(gridPointSorting) {
  assert(_gridFactory);
  _geometry->addSensitiveObject(this->_self);
}

std::unique_ptr<AtomCenteredGrid> AtomCenteredGridController::getAtomGrid() {
  return _gridFactory->produce(_geometry);
}

void AtomCenteredGridController::produceGrid() {
  auto newGrid = _gridFactory->produce(_geometry);
  _grid = std::move(newGrid);
  if (_gridPointSorting)
    _grid->sort();
}

void AtomCenteredGridController::notify() {
  _grid = nullptr;
  notifyObjects();
}

} /* namespace Serenity */
