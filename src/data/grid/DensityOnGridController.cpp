/**
 * @file   DensityOnGridController.cpp
 * @author Thomas Dresselhaus
 *
 * @date   22. November 2014, 10:55
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
#include "data/grid/DensityOnGridController.h"
/* Include Serenity Internal Headers */
#include "grid/GridController.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DensityOnGridController<SCFMode>::DensityOnGridController(std::shared_ptr<GridController> gridController,
                                                          const unsigned int highestDerivative)
  : _gridController(gridController),
    _highestDerivative(highestDerivative),
    _nGridPoints(gridController->getNGridPoints()),
    _upToDate(false) {
  assert(_gridController);
  _gridController->addSensitiveObject(this->_self);
}

template<Options::SCF_MODES SCFMode>
unsigned int DensityOnGridController<SCFMode>::getNGridPoints() const {
  return _gridController->getNGridPoints();
}

template class DensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class DensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity
