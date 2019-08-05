/**
 * @file   BasisFunctionOnGridControllerFactory.cpp
 *
 * @date   May 8, 2014
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
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "input/FunctionalClassResolver.h"
#include "grid/GridControllerFactory.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <cassert>


namespace Serenity {
using namespace std;

/*
 * Initialize static pointer for the singleton instance.
 */
std::unique_ptr<BasisFunctionOnGridControllerFactory>
    BasisFunctionOnGridControllerFactory::_instance;

shared_ptr<BasisFunctionOnGridController> BasisFunctionOnGridControllerFactory::produce(
    unsigned int maxBlockSize,
    double radialThreshold,
    unsigned int highestDerivative,
    shared_ptr<BasisController> basisController,
    shared_ptr<GridController> gridController) {
  if (!_instance) _instance.reset(new BasisFunctionOnGridControllerFactory());
  return _instance->getOrProduce(
      basisController, gridController, maxBlockSize, radialThreshold, highestDerivative);
}

shared_ptr<BasisFunctionOnGridController> BasisFunctionOnGridControllerFactory::produce(
    const Settings& settings,
    shared_ptr<BasisController> basisController,
    shared_ptr<GridController> gridController,
    unsigned int highestDerivative) {
  return produce(
      settings.grid.blocksize,
      settings.grid.basFuncRadialThreshold,
      highestDerivative,
      basisController,
      gridController);
}

unique_ptr<BasisFunctionOnGridController> BasisFunctionOnGridControllerFactory::produceNew(
      const std::shared_ptr<BasisController> basisController,
      const std::shared_ptr<GridController> gridController,
      const unsigned int maxBlockSize,
      const double radialThreshold,
      const unsigned int highestDerivative) {
  return unique_ptr<BasisFunctionOnGridController>(new BasisFunctionOnGridController(
      basisController,
      gridController,
      maxBlockSize,
      radialThreshold,
      highestDerivative));
}

} /* namespace Serenity */
