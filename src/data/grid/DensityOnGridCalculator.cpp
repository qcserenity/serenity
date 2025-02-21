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

template<Options::SCF_MODES SCFMode>
DensityOnGridCalculator<SCFMode>::DensityOnGridCalculator(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                          const double blockAverageThreshold)
  : _basisFunctionOnGridController(basisFunctionOnGridController), _blockAverageThreshold(blockAverageThreshold) {
  assert(_basisFunctionOnGridController);
}

template<Options::SCF_MODES SCFMode>
void DensityOnGridCalculator<SCFMode>::calcDensityOnGrid(const DensityMatrix<SCFMode>& densityMatrix,
                                                         DensityOnGrid<SCFMode>& densityOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));
  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES SCFMode>
DensityOnGrid<SCFMode> DensityOnGridCalculator<SCFMode>::calcDensityOnGrid(const DensityMatrix<SCFMode>& densityMatrix) {
  DensityOnGrid<SCFMode> densOnGrid(_basisFunctionOnGridController->getGridController());
  calcDensityOnGrid(densityMatrix, densOnGrid);
  return densOnGrid;
}

template<Options::SCF_MODES SCFMode>
void DensityOnGridCalculator<SCFMode>::calcDensityAndGradientOnGrid(const DensityMatrix<SCFMode>& densityMatrix,
                                                                    DensityOnGrid<SCFMode>& densityOnGrid,
                                                                    Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));

  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                              *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES SCFMode>
DensityOnGrid<SCFMode>
DensityOnGridCalculator<SCFMode>::calcDensityAndGradientOnGrid(const DensityMatrix<SCFMode>& densityMatrix,
                                                               Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid) {
  DensityOnGrid<SCFMode> densOnGrid(_basisFunctionOnGridController->getGridController());
  calcDensityAndGradientOnGrid(densityMatrix, densOnGrid, densityGradientOnGrid);
  return densOnGrid;
}

template<Options::SCF_MODES SCFMode>
void DensityOnGridCalculator<SCFMode>::calcDensityAndDerivativesOnGrid(const DensityMatrix<SCFMode>& densityMatrix,
                                                                       DensityOnGrid<SCFMode>& densityOnGrid,
                                                                       Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid,
                                                                       Hessian<DensityOnGrid<SCFMode>>& densityHessianOnGrid) {
  Timings::takeTime("Tech. -  Density On Grid Eval.");
  assert(isDefinedInSameBasis(densityMatrix, *_basisFunctionOnGridController));
  assert(isDefinedOnSameGrid(densityOnGrid, *_basisFunctionOnGridController));

  _nonNegligible = MatrixOperatorToGridTransformer::transform(densityMatrix, densityOnGrid, densityGradientOnGrid,
                                                              densityHessianOnGrid, *_basisFunctionOnGridController);
  Timings::timeTaken("Tech. -  Density On Grid Eval.");
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<GridController> DensityOnGridCalculator<SCFMode>::getGridController() const {
  return _basisFunctionOnGridController->getGridController();
}

template class DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>;
template class DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
