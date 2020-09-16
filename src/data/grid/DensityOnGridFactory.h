/**
 * @file   DensityOnGridFactory.h
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

#ifndef DATA_GRID_DENSITYONGRIDFACTORY_H_
#define DATA_GRID_DENSITYONGRIDFACTORY_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "misc/RememberingFactory.h"

namespace Serenity {

class Geometry;
class GridController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
struct Settings;

/**
 * @class DensityOnGridFactory DensityOnGridFactory.h
 * @brief A factory for DensityOnGridController, in order to avoid multiple constructions
 *          and recalculation of DensityOnGrid type data.
 */
template<Options::SCF_MODES SCFMode>
class DensityOnGridFactory
  : public RememberingFactory<DensityOnGridController<SCFMode>, const std::shared_ptr<DensityMatrixController<SCFMode>>,
                              const std::shared_ptr<GridController>, const unsigned int, const unsigned int, const double, const double> {
 private:
  /**
   * @brief Private default constructor - Singleton.
   */
  DensityOnGridFactory() = default;

 public:
  /// @brief Default destructor.
  virtual ~DensityOnGridFactory() = default;

  /**
   * @brief For a detailed description see DensityOnGridController.
   * @param density
   * @param radialThreshold
   * @param gridController
   * @param highestDerivative
   * @param settings
   * @return Returns a DensityOnGridController.
   */
  static std::shared_ptr<DensityOnGridController<SCFMode>>
  produce(const std::shared_ptr<DensityMatrixController<SCFMode>> density, const std::shared_ptr<GridController> gridController,
          const unsigned int highestDerivative, const Settings& settings);

 private:
  std::unique_ptr<DensityOnGridController<SCFMode>>
  produceNew(const std::shared_ptr<DensityMatrixController<SCFMode>> density,
             const std::shared_ptr<GridController> gridController, const unsigned int highestDerivative,
             const unsigned int blocksize, const double basFuncRadialThreshold, const double blockAveThreshold) override final;

  static std::unique_ptr<DensityOnGridFactory<SCFMode>> _instance;
};

} /* namespace Serenity */

#endif /* DATA_GRID_DENSITYONGRIDFACTORY_H_ */
