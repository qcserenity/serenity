/**
 * @file   DensityMatrixDensityOnGridController.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 * 
 * @date   30. Dezember 2014, 16:05
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
#include "data/grid/DensityMatrixDensityOnGridController.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "math/FloatMaths.h"
#include "io/FormattedOutput.h"
#include "grid/GridController.h"
#include "io/IOOptions.h"
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <cassert>
#include <stdexcept>
#include <string>


namespace Serenity {

template<Options::SCF_MODES T>
DensityMatrixDensityOnGridController<T>::DensityMatrixDensityOnGridController(
      std::shared_ptr<DensityOnGridCalculator<T> > densOnGridCalculator,
      const std::shared_ptr<DensityMatrixController<T> > densityMatrixController) :
    DensityMatrixDensityOnGridController<T>(
      densOnGridCalculator,
      densityMatrixController,
      densOnGridCalculator->getBasisFunctionOnGridController()->getHighestDerivative()) {
}

template<Options::SCF_MODES T>
DensityMatrixDensityOnGridController<T>::DensityMatrixDensityOnGridController(
      std::shared_ptr<DensityOnGridCalculator<T> > densOnGridCalculator,
      const std::shared_ptr<DensityMatrixController<T> > densityMatrixController,
      const unsigned int highestDerivative) :
    DensityOnGridController<T>(densOnGridCalculator->getGridController(), highestDerivative),
    _densOnGridCalculator(densOnGridCalculator),
    _densityMatrixController(densityMatrixController) {
  assert(_densOnGridCalculator);
  _densityMatrixController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<T> >::_self);
  this->_densityOnGrid.reset(new DensityOnGrid<T>(this->getGridController()));
  if (this->_highestDerivative >= 1)
    this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<T> >(this->getGridController());
  if (this->_highestDerivative >= 2)
    this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<T> >(this->getGridController());
  assert(this->_highestDerivative <=
      _densOnGridCalculator->getBasisFunctionOnGridController()->getHighestDerivative());
}

template<Options::SCF_MODES T>
const DensityOnGrid<T>& DensityMatrixDensityOnGridController<T>::getDensityOnGrid() {
  if (!this->_upToDate) calculateData();
  assert(this->_densityOnGrid->isValid());
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES T>
const Gradient<DensityOnGrid<T> >&
    DensityMatrixDensityOnGridController<T>::getDensityGradientOnGrid() {
//  assert(this->_highestDerivative >= 1);
  if (this->_highestDerivative < 1) this->setHighestDerivative(1);
  if (!this->_upToDate) calculateData();
  for (const auto& component : *this->_densityGradientOnGrid){
    if (!component.isValid())
      throw SerenityError("A component of the Density stored on the grid is invalid.");
  }
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES T>
const Hessian<DensityOnGrid<T> >&
    DensityMatrixDensityOnGridController<T>::getDensityHessianOnGrid() {
  if (this->_highestDerivative < 2) this->setHighestDerivative(2);
  if (!this->_upToDate) calculateData();
  for (const auto& component : *this->_densityHessianOnGrid){
    if (!component.isValid())
      throw SerenityError("A component of the Density stored on the grid is invalid.");
  }
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES T>void DensityMatrixDensityOnGridController<T>::notify() {
  this->_upToDate = false;
  this->notifyObjects();
}

template<Options::SCF_MODES T>void DensityMatrixDensityOnGridController<T>::setHighestDerivative(
    unsigned int newHighestDerivative) {
  assert(newHighestDerivative <= 2);
  if (newHighestDerivative > this->_highestDerivative) {
    if (this->_highestDerivative < 1 && newHighestDerivative >= 1)
      this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<T> >(this->getGridController());
    if (this->_highestDerivative < 2 && newHighestDerivative >= 2)
      this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<T> >(this->getGridController());
    if (this->_upToDate)
      std::cout << "Warning! A new highest derivative is set in DensityMatrixDensityOnGridController "
                   "causing that data is thrown away and partly recalculated. This is inefficient "
                   "and should be fixed in the code!" << std::endl;
    // TODO actually not needed like this; data for other objects may still be valid
    this->notify();
  } else if (newHighestDerivative < this->_highestDerivative) {
    if (newHighestDerivative < 2 && this->_highestDerivative >= 2)
      this->_densityHessianOnGrid.reset();
    if (newHighestDerivative < 1 && this->_highestDerivative >= 1)
      this->_densityGradientOnGrid.reset();
    this->notify();
  }
  this->_highestDerivative = newHighestDerivative;
}

template<Options::SCF_MODES T>void DensityMatrixDensityOnGridController<T>::calculateData() {
  auto densityMatrix = _densityMatrixController->getDensityMatrix();

  switch(this->_highestDerivative) {
    case 0:
      _densOnGridCalculator->calcDensityOnGrid(densityMatrix, *this->_densityOnGrid);
      break;
    case 1:
      _densOnGridCalculator->calcDensityAndGradientOnGrid(
          densityMatrix, *this->_densityOnGrid, *this->_densityGradientOnGrid);
      break;
    case 2:
      _densOnGridCalculator->calcDensityAndDerivativesOnGrid(
          densityMatrix,
          *this->_densityOnGrid,
          *this->_densityGradientOnGrid,
          *this->_densityHessianOnGrid);
      break;
    default:
      throw SerenityError("Derivative of density on grid higher than 2 was requested. Not implemented.");
  }
  /*
   * Check for the grid accuracy
   */
  if (iOOptions.gridAccuracyCheck) {
    const auto& weights = this->getGridController()->getWeights();
    SpinPolarizedData<T, double> nElectrons(0.0);
    const auto& dens = *this->_densityOnGrid;
    for_spin(nElectrons, dens) {
      nElectrons_spin += dens_spin.dot(weights);
      print((std::string)"nElectrons from integration over grid: " + nElectrons_spin);
    };
  }
  this->_upToDate = true;
}

template<Options::SCF_MODES T>
void DensityMatrixDensityOnGridController<T>::setDensityOnGrid(
    std::unique_ptr<DensityOnGrid<T>> densityOnGrid){
  this->_densityOnGrid.reset(densityOnGrid.release());
  this->_upToDate=true;
};

template<Options::SCF_MODES T>
void DensityMatrixDensityOnGridController<T>::setDensityGradientOnGrid(
    std::unique_ptr<Gradient<DensityOnGrid<T> >> densityGradientOnGrid){
  this->_densityGradientOnGrid.reset(densityGradientOnGrid.release());
  this->_upToDate=true;
};

template<Options::SCF_MODES T>
void DensityMatrixDensityOnGridController<T>::setDensityHessianOnGrid(
    std::unique_ptr<Hessian<DensityOnGrid<T> >> densityHessianOnGrid){
  this->_densityHessianOnGrid.reset(densityHessianOnGrid.release());
  this->_upToDate=true;
};

template class DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
