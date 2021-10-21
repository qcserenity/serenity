/**
 * @file   DensityOnGridFactory.cpp
 *
 * @date   Oct 19, 2017
 * @author Jan Unsleber
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
#include "data/grid/DensityOnGridFactory.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrixController.h"
#include "grid/GridController.h"
#include "settings/Settings.h"

namespace Serenity {

/*
 * Initialize static pointer for the singleton instance.
 */
template<Options::SCF_MODES SCFMode>
std::unique_ptr<DensityOnGridFactory<SCFMode>> DensityOnGridFactory<SCFMode>::_instance;

template<Options::SCF_MODES SCFMode>
std::shared_ptr<DensityOnGridController<SCFMode>>
DensityOnGridFactory<SCFMode>::produce(const std::shared_ptr<DensityMatrixController<SCFMode>> density,
                                       const std::shared_ptr<GridController> gridController,
                                       const unsigned int highestDerivative, const Settings& settings) {
  if (!_instance)
    _instance.reset(new DensityOnGridFactory<SCFMode>());
  return _instance->getOrProduce(density, gridController, highestDerivative, settings.grid.blocksize,
                                 settings.grid.basFuncRadialThreshold, settings.grid.blockAveThreshold);
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<DensityOnGridController<SCFMode>>
DensityOnGridFactory<SCFMode>::produceNew(const std::shared_ptr<DensityMatrixController<SCFMode>> density,
                                          const std::shared_ptr<GridController> gridController,
                                          const unsigned int highestDerivative, const unsigned int blocksize,
                                          const double basFuncRadialThreshold, const double blockAveThreshold) {
  auto basisFunctionOnGridController(BasisFunctionOnGridControllerFactory::produce(
      blocksize, basFuncRadialThreshold, highestDerivative, density->getDensityMatrix().getBasisController(), gridController));
  auto densityOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<SCFMode>>(basisFunctionOnGridController, blockAveThreshold);
  auto densOnGrid = std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, density);

  return std::make_unique<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, density);
}

template class DensityOnGridFactory<Options::SCF_MODES::RESTRICTED>;
template class DensityOnGridFactory<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
