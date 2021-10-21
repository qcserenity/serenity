/**
 * @file   DensityOnGridCalculator.cpp
 *
 * @date   Mar 15, 2014
 * @author Dennis Barton, severely modified by Thomas Dresselhaus
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
#include "data/grid/DensityOnGridCalculator.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/MatrixOperatorToGridTransformer.h"
#include "data/matrices/DensityMatrix.h"
#include "misc/Timing.h"

namespace Serenity {

template<Options::SCF_MODES T>
DensityOnGridCalculator<T>::DensityOnGridCalculator(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                    const double blockAverageThreshold)
  : _basisFunctionOnGridController(basisFunctionOnGridController), _blockAverageThreshold(blockAverageThreshold) {
  assert(_basisFunctionOnGridController);
}

template<Options::SCF_MODES T>
void DensityOnGridCalculator<T>::calcDensityOnGrid(const DensityMatrix<T>& densityMatrix, DensityOnGrid<T>& densityOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));
  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES T>
DensityOnGrid<T> DensityOnGridCalculator<T>::calcDensityOnGrid(const DensityMatrix<T>& densityMatrix) {
  DensityOnGrid<T> densOnGrid(_basisFunctionOnGridController->getGridController());
  calcDensityOnGrid(densityMatrix, densOnGrid);
  return densOnGrid;
}

template<Options::SCF_MODES T>
void DensityOnGridCalculator<T>::calcDensityAndGradientOnGrid(const DensityMatrix<T>& densityMatrix,
                                                              DensityOnGrid<T>& densityOnGrid,
                                                              Gradient<DensityOnGrid<T>>& densityGradientOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));

  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                              *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES T>
DensityOnGrid<T> DensityOnGridCalculator<T>::calcDensityAndGradientOnGrid(const DensityMatrix<T>& densityMatrix,
                                                                          Gradient<DensityOnGrid<T>>& densityGradientOnGrid) {
  DensityOnGrid<T> densOnGrid(_basisFunctionOnGridController->getGridController());
  calcDensityAndGradientOnGrid(densityMatrix, densOnGrid, densityGradientOnGrid);
  return densOnGrid;
}

template<Options::SCF_MODES T>
void DensityOnGridCalculator<T>::calcDensityAndDerivativesOnGrid(const DensityMatrix<T>& densityMatrix,
                                                                 DensityOnGrid<T>& densityOnGrid,
                                                                 Gradient<DensityOnGrid<T>>& densityGradientOnGrid,
                                                                 Hessian<DensityOnGrid<T>>& densityHessianOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));

  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                              densityHessianOnGrid, *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES T>
std::shared_ptr<GridController> DensityOnGridCalculator<T>::getGridController() const {
  return _basisFunctionOnGridController->getGridController();
}

template class DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>;
template class DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
