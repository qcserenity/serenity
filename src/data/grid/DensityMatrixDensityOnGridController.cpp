/**
 * @file   DensityMatrixDensityOnGridController.cpp
 * @author Thomas Dresselhaus
 *
 * @date   Dec 30, 2014
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
#include "data/grid/DensityMatrixDensityOnGridController.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
/* Include Std and External Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DensityMatrixDensityOnGridController<SCFMode>::DensityMatrixDensityOnGridController(
    std::shared_ptr<DensityOnGridCalculator<SCFMode>> densOnGridCalculator,
    const std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController)
  : DensityMatrixDensityOnGridController<SCFMode>(
        densOnGridCalculator, densityMatrixController,
        densOnGridCalculator->getBasisFunctionOnGridController()->getHighestDerivative()) {
}

template<Options::SCF_MODES SCFMode>
DensityMatrixDensityOnGridController<SCFMode>::DensityMatrixDensityOnGridController(
    std::shared_ptr<DensityOnGridCalculator<SCFMode>> densOnGridCalculator,
    const std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController, const unsigned int highestDerivative)
  : DensityOnGridController<SCFMode>(densOnGridCalculator->getGridController(), highestDerivative),
    _densOnGridCalculator(densOnGridCalculator),
    _densityMatrixController(densityMatrixController) {
  assert(_densOnGridCalculator);
  _densityMatrixController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  this->_densityOnGrid.reset(new DensityOnGrid<SCFMode>(this->getGridController()));
  if (this->_highestDerivative >= 1) {
    this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(this->getGridController());
  }
  if (this->_highestDerivative >= 2) {
    this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<SCFMode>>(this->getGridController());
  }
  assert(this->_highestDerivative <= _densOnGridCalculator->getBasisFunctionOnGridController()->getHighestDerivative());
}

template<Options::SCF_MODES SCFMode>
const DensityOnGrid<SCFMode>& DensityMatrixDensityOnGridController<SCFMode>::getDensityOnGrid() {
  if (!this->_upToDate) {
    this->calculateData();
  }
  assert(this->_densityOnGrid->isValid());
  return *this->_densityOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Gradient<DensityOnGrid<SCFMode>>& DensityMatrixDensityOnGridController<SCFMode>::getDensityGradientOnGrid() {
  if (this->_highestDerivative < 1) {
    this->setHighestDerivative(1);
  }
  if (!this->_upToDate) {
    this->calculateData();
  }
  for (const auto& component : *this->_densityGradientOnGrid) {
    if (!component.isValid()) {
      throw SerenityError("A component of the Density stored on the grid is invalid.");
    }
  }
  return *this->_densityGradientOnGrid;
}

template<Options::SCF_MODES SCFMode>
const Hessian<DensityOnGrid<SCFMode>>& DensityMatrixDensityOnGridController<SCFMode>::getDensityHessianOnGrid() {
  if (this->_highestDerivative < 2) {
    this->setHighestDerivative(2);
  }
  if (!this->_upToDate) {
    this->calculateData();
  }
  for (const auto& component : *this->_densityHessianOnGrid) {
    if (!component.isValid()) {
      throw SerenityError("A component of the Density stored on the grid is invalid.");
    }
  }
  return *this->_densityHessianOnGrid;
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::notify() {
  this->_upToDate = false;
  this->notifyObjects();
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::setHighestDerivative(unsigned int newHighestDerivative) {
  assert(newHighestDerivative <= 2);
  if (newHighestDerivative > this->_highestDerivative) {
    if (this->_highestDerivative < 1 && newHighestDerivative >= 1)
      this->_densityGradientOnGrid = makeGradientPtr<DensityOnGrid<SCFMode>>(this->getGridController());
    if (this->_highestDerivative < 2 && newHighestDerivative >= 2)
      this->_densityHessianOnGrid = makeHessianPtr<DensityOnGrid<SCFMode>>(this->getGridController());
    if (this->_upToDate)
      std::cout << "Warning! A new highest derivative is set in DensityMatrixDensityOnGridController "
                   "causing that data is thrown away and partly recalculated. This is inefficient "
                   "and should be fixed in the code!"
                << std::endl;
    // TODO actually not needed like this; data for other objects may still be valid
    this->notify();
  }
  else if (newHighestDerivative < this->_highestDerivative) {
    if (newHighestDerivative < 2 && this->_highestDerivative >= 2) {
      this->_densityHessianOnGrid.reset();
    }
    if (newHighestDerivative < 1 && this->_highestDerivative >= 1) {
      this->_densityGradientOnGrid.reset();
    }
    this->notify();
  }
  this->_highestDerivative = newHighestDerivative;
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::calculateData() {
  auto densityMatrix = _densityMatrixController->getDensityMatrix();

  if (this->_highestDerivative == 0) {
    _densOnGridCalculator->calcDensityOnGrid(densityMatrix, *this->_densityOnGrid);
  }
  else if (this->_highestDerivative == 1) {
    _densOnGridCalculator->calcDensityAndGradientOnGrid(densityMatrix, *this->_densityOnGrid, *this->_densityGradientOnGrid);
  }
  else if (this->_highestDerivative == 2) {
    _densOnGridCalculator->calcDensityAndDerivativesOnGrid(densityMatrix, *this->_densityOnGrid,
                                                           *this->_densityGradientOnGrid, *this->_densityHessianOnGrid);
  }
  else {
    throw SerenityError("Derivative of density on grid higher than 2 was requested. Not implemented.");
  }

  // Check for the grid accuracy
  if (iOOptions.gridAccuracyCheck) {
    const auto& weights = this->getGridController()->getWeights();
    SpinPolarizedData<SCFMode, double> nElectrons(0.0);
    const auto& dens = *this->_densityOnGrid;
    for_spin(nElectrons, dens) {
      nElectrons_spin += dens_spin.dot(weights);
      print((std::string) "nElectrons from integration over grid: " + nElectrons_spin);
    };
  }
  this->_upToDate = true;
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::setDensityOnGrid(std::unique_ptr<DensityOnGrid<SCFMode>> densityOnGrid) {
  this->_densityOnGrid = std::move(densityOnGrid);
  this->_upToDate = true;
};

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::setDensityGradientOnGrid(
    std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> densityGradientOnGrid) {
  this->_densityGradientOnGrid = std::move(densityGradientOnGrid);
  this->_upToDate = true;
};

template<Options::SCF_MODES SCFMode>
void DensityMatrixDensityOnGridController<SCFMode>::setDensityHessianOnGrid(
    std::unique_ptr<Hessian<DensityOnGrid<SCFMode>>> densityHessianOnGrid) {
  this->_densityHessianOnGrid = std::move(densityHessianOnGrid);
  this->_upToDate = true;
};

template class DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>;
template class DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
