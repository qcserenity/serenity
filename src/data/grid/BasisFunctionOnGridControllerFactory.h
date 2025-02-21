/**
 * @file   BasisFunctionOnGridControllerFactory.h
 *
 * @date   May 8, 2014
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
#ifndef BASISFUNCTIONONGRIDCONTROLLERFACTORY_H_
#define BASISFUNCTIONONGRIDCONTROLLERFACTORY_H_
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "misc/RememberingFactory.h"

namespace Serenity {
/* Forward declarations */
struct Settings;
/**
 * @class BasisFunctionOnGridControllerFactory BasisFunctionOnGridControllerFactory.h
 * @brief Produces instances of BasisFunctionOnGridController
 */
class BasisFunctionOnGridControllerFactory
  : public RememberingFactory<BasisFunctionOnGridController, const std::shared_ptr<BasisController>,
                              const std::shared_ptr<GridController>, const unsigned int, const double, const unsigned int> {
 private:
  /**
   * private default constructor - Singleton
   */
  BasisFunctionOnGridControllerFactory() = default;

 public:
  virtual ~BasisFunctionOnGridControllerFactory() = default;
  /**
   * @see BasisFunctionOnGridController
   * @param maxBlockSize The maximum number of grid points per block.
   * @param radialThreshold Basis functions with a radial value below this threshold
   *        will be skipped.
   * @param highestDerivative The highest derivate d^n phi/(d^n r) to be computed.
   * @param basisController Provides the basis functions.
   * @param gridController Provides the used grid.
   * @return A BasisFunctionOnGridController.
   */
  static std::shared_ptr<BasisFunctionOnGridController>
  produce(unsigned int maxBlockSize, double radialThreshold, unsigned int highestDerivative,
          std::shared_ptr<BasisController> basisController, std::shared_ptr<GridController> gridController);
  /**
   * @param settings Contains values for maxBlockSize and radialThreshold
   *                 from user input.
   * @param basisController Provides the basis functions.
   * @param gridController Provides the used grid.
   * @param highestDerivatve The highest derivate d^n phi/(d^n r) to be computed.
   * @return A BasisFunctionOnGridController.
   */
  static std::shared_ptr<BasisFunctionOnGridController>
  produce(const Settings& settings, std::shared_ptr<BasisController> basisController,
          std::shared_ptr<GridController> gridController, unsigned int highestDerivative = 1);

 private:
  std::unique_ptr<BasisFunctionOnGridController> produceNew(const std::shared_ptr<BasisController> basis,
                                                            const std::shared_ptr<GridController> gridController,
                                                            const unsigned int maxBlockSize, const double radialThreshold,
                                                            const unsigned int highestDerivative) override final;

  static std::unique_ptr<BasisFunctionOnGridControllerFactory> _instance;
};

} /* namespace Serenity */

#endif /* BASISFUNCTIONONGRIDCONTROLLERFACTORY_H_ */
