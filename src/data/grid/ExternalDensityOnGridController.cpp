/**
 * @file ExternalDensityOnGridController.cpp
 * 
 * @date Dezember 30, 2014
 * @author Thomas Dresselhaus
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
/* Include Class Header*/
#include "data/grid/ExternalDensityOnGridController.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {

template<Options::SCF_MODES T>
ExternalDensityOnGridController<T>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<T> >& densityOnGrid) :
      DensityOnGridController<T>(densityOnGrid->getGridController(), 0) {
  this->_densityOnGrid.swap(densityOnGrid);
}

template<Options::SCF_MODES T>
ExternalDensityOnGridController<T>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<T> >& densityOnGrid,
    std::unique_ptr<Gradient<DensityOnGrid<T> > >& densityGradientOnGrid) :
      DensityOnGridController<T>(densityOnGrid->getGridController(), 1) {
  this->_densityOnGrid.swap(densityOnGrid);
  this->_densityGradientOnGrid.swap(densityGradientOnGrid);
  for (const auto& component : *this->_densityGradientOnGrid){
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
}

template<Options::SCF_MODES T>
ExternalDensityOnGridController<T>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<T> >& densityOnGrid,
    std::unique_ptr<Gradient<DensityOnGrid<T> > >& densityGradientOnGrid,
    std::unique_ptr<Hessian<DensityOnGrid<T> > >& densityHessianOnGrid) :
      DensityOnGridController<T>(densityOnGrid->getGridController(), 2) {
  this->_densityOnGrid.swap(densityOnGrid);
  this->_densityGradientOnGrid.swap(densityGradientOnGrid);
  this->_densityHessianOnGrid.swap(densityHessianOnGrid);
  for (const auto& component : *this->_densityGradientOnGrid){
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
  for (const auto& component : *this->_densityHessianOnGrid){
    if (not isDefinedOnSameGrid(component, *this->_densityOnGrid))
      throw SerenityError("DensityOnGridController: data is not defined on the same grind");
  }
}

template<Options::SCF_MODES T>
ExternalDensityOnGridController<T>::ExternalDensityOnGridController(
    std::unique_ptr<DensityOnGrid<T> >& densityOnGrid,
    std::unique_ptr<Gradient<DensityOnGrid<T> > >& densityGradientOnGrid,
    std::unique_ptr<Hessian<DensityOnGrid<T> > >& densityHessianOnGrid,
    std::function<void(DensityOnGrid<T>& densityOnGrid,
                       Gradient<DensityOnGrid<T> >& densityGradientOnGrid,
                       Hessian<DensityOnGrid<T> >& densityHessianOnGrid)> externalUpdateFunction) :
  ExternalDensityOnGridController<T>(densityOnGrid, densityGradientOnGrid, densityHessianOnGrid) {
  // The constructor delegation somehow makes initializing _externalUpdateFunction directly impossible...
    _externalUpdateFunction = externalUpdateFunction;
}

template<Options::SCF_MODES T>
ExternalDensityOnGridController<T>::~ExternalDensityOnGridController() {
}

template<Options::SCF_MODES T>
const DensityOnGrid<T>& ExternalDensityOnGridController<T>::getDensityOnGrid() {
  assert(this->_densityOnGrid->isValid());
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES T>
const Gradient<DensityOnGrid<T> >& ExternalDensityOnGridController<T>::getDensityGradientOnGrid() {
  for (const auto& component : *this->_densityGradientOnGrid){
    if (!component.isValid())
      throw SerenityError("A component of the Density stored on the grid is invalid.");
  }
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES T>
const Hessian<DensityOnGrid<T> >& ExternalDensityOnGridController<T>::getDensityHessianOnGrid() {
  for (const auto& component : *this->_densityHessianOnGrid){
    if (!component.isValid())
      throw SerenityError("A component of the Density stored on the grid is invalid.");
  }
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES T>
void ExternalDensityOnGridController<T>::notify() {
  if (_externalUpdateFunction) {
    _externalUpdateFunction(
        *this->_densityOnGrid, *this->_densityGradientOnGrid, *this->_densityHessianOnGrid);
    this->notifyObjects();
  }
}

template class ExternalDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class ExternalDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
