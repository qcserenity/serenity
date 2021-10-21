/**
 * @file DensityAdder.cpp
 *
 * @date Sep 5, 2016
 * @author Michael Boeckers
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
#include "data/grid/DensityAdder.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void DensityAdder<SCFMode>::add(DensityOnGrid<SCFMode>& densityToBeAddedTo, const DensityMatrix<SCFMode>& densityMatrix,
                                std::shared_ptr<BasisController> basisController, const double blockAverageThreshold) {
  // get grid controller of grid on which the result density is defined
  auto gridController = densityToBeAddedTo.getGridController();
  // get basisFunctionOnGridController for result grid
  auto basisFunctionOnGridController =
      std::make_shared<BasisFunctionOnGridController>(basisController, gridController, 128, 1e-9, 0);
  // get density on grid calculator to calculate density on result grid
  DensityOnGridCalculator<SCFMode> densityOnGridCalculator(basisFunctionOnGridController, blockAverageThreshold);
  // calculate density on result grid
  DensityOnGrid<SCFMode> densityToBeAdded(gridController);
  densityOnGridCalculator.calcDensityOnGrid(densityMatrix, densityToBeAdded);
  // Add densities
  densityToBeAddedTo += densityToBeAdded;
}

template class DensityAdder<Options::SCF_MODES::RESTRICTED>;
template class DensityAdder<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
