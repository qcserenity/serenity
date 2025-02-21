/**
 * @file ExternalDensityOnGridController.cpp
 *
 * @date Dezember 30, 2014
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
#include "data/grid/ExternalDensityOnGridController.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ExternalDensityOnGridController<SCFMode>::ExternalDensityOnGridController(std::unique_ptr<DensityOnGrid<SCFMode>>& densityOnGrid)
  : DensityOnGridController<SCFMode>(densityOnGrid->getGridController(), 0) {
  this->_densityOnGrid = std::move(densityOnGrid);
}

template<Options::SCF_MODES SCFMode>
ExternalDensityOnGridController<SCFMode>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<SCFMode>>& densityOnGrid, std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>>& densityGradientOnGrid)
  : DensityOnGridController<SCFMode>(densityOnGrid->getGridController(), 1) {
  this->_densityOnGrid = std::move(densityOnGrid);
  this->_densityGradientOnGrid = std::move(densityGradientOnGrid);
  for (const auto& component : *this->_densityGradientOnGrid) {
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
}

template<Options::SCF_MODES SCFMode>
ExternalDensityOnGridController<SCFMode>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<SCFMode>>& densityOnGrid,
    std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>>& densityGradientOnGrid,
    std::unique_ptr<Hessian<DensityOnGrid<SCFMode>>>& densityHessianOnGrid)
  : DensityOnGridController<SCFMode>(densityOnGrid->getGridController(), 2) {
  this->_densityOnGrid = std::move(densityOnGrid);
  this->_densityGradientOnGrid = std::move(densityGradientOnGrid);
  this->_densityHessianOnGrid = std::move(densityHessianOnGrid);
  for (const auto& component : *this->_densityGradientOnGrid) {
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
  for (const auto& component : *this->_densityHessianOnGrid) {
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
}

template<Options::SCF_MODES SCFMode>
ExternalDensityOnGridController<SCFMode>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<SCFMode>>& densityOnGrid,
    std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>>& densityGradientOnGrid,
    std::unique_ptr<Hessian<DensityOnGrid<SCFMode>>>& densityHessianOnGrid,
    std::function<void(DensityOnGrid<SCFMode>& densityOnGrid, Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid,
                       Hessian<DensityOnGrid<SCFMode>>& densityHessianOnGrid)>
        externalUpdateFunction)
  : ExternalDensityOnGridController<SCFMode>(densityOnGrid, densityGradientOnGrid, densityHessianOnGrid) {
  // The constructor delegation somehow makes initializing _externalUpdateFunction directly impossible...
  _externalUpdateFunction = externalUpdateFunction;
}

template<Options::SCF_MODES SCFMode>
ExternalDensityOnGridController<SCFMode>::~ExternalDensityOnGridController() {
}

template<Options::SCF_MODES SCFMode>
const DensityOnGrid<SCFMode>& ExternalDensityOnGridController<SCFMode>::getDensityOnGrid() {
  if (!this->_densityOnGrid->isValid()) {
    throw SerenityError("A component of the density stored on the grid is invalid.");
  }
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Gradient<DensityOnGrid<SCFMode>>& ExternalDensityOnGridController<SCFMode>::getDensityGradientOnGrid() {
  for (const auto& component : *this->_densityGradientOnGrid) {
    if (!component.isValid())
      throw SerenityError("A component of the density gradient stored on the grid is invalid.");
  }
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Hessian<DensityOnGrid<SCFMode>>& ExternalDensityOnGridController<SCFMode>::getDensityHessianOnGrid() {
  for (const auto& component : *this->_densityHessianOnGrid) {
    if (!component.isValid())
      throw SerenityError("A component of the density hessian stored on the grid is invalid.");
  }
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES SCFMode>
void ExternalDensityOnGridController<SCFMode>::notify() {
  if (_externalUpdateFunction) {
    _externalUpdateFunction(*this->_densityOnGrid, *this->_densityGradientOnGrid, *this->_densityHessianOnGrid);
    this->notifyObjects();
  }
}

template class ExternalDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class ExternalDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
