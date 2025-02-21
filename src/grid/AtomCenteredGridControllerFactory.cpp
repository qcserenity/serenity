/**
 * @file   AtomCenteredGridControllerFactory.cpp
 *
 * @date   Apr 29, 2014
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
/* Include Class Header*/
#include "grid/AtomCenteredGridControllerFactory.h"
/* Include Serenity Internal Headers */
#include "grid/construction/GridFactory.h"
#include "settings/Settings.h"

namespace Serenity {

std::shared_ptr<AtomCenteredGridController>
AtomCenteredGridControllerFactory::produce(std::shared_ptr<const Geometry> geometry, Options::GRID_TYPES gridType,
                                           unsigned int gridSmoothing, Options::RADIAL_GRID_TYPES radialGridType,
                                           Options::SPHERICAL_GRID_TYPES sphericalGridType, unsigned int gridAccuracy,
                                           double gridWeightThreshold, const bool gridPointSorting) {
  if (!_instance)
    _instance.reset(new AtomCenteredGridControllerFactory);
  return _instance->getOrProduce(geometry, gridType, gridSmoothing, radialGridType, sphericalGridType, gridAccuracy,
                                 gridWeightThreshold, gridPointSorting);
}

std::shared_ptr<AtomCenteredGridController> AtomCenteredGridControllerFactory::produce(std::shared_ptr<const Geometry> geometry,
                                                                                       const GRID& settings,
                                                                                       Options::GRID_PURPOSES purpose) {
  return produce(geometry, settings.gridType, settings.smoothing, settings.radialGridType, settings.sphericalGridType,
                 (purpose == Options::GRID_PURPOSES::SMALL) ? settings.smallGridAccuracy : settings.accuracy,
                 settings.weightThreshold, settings.gridPointSorting);
}

std::unique_ptr<AtomCenteredGridController>
AtomCenteredGridControllerFactory::produceNew(const std::shared_ptr<const Geometry> geometry, Options::GRID_TYPES gridType,
                                              unsigned int gridSmoothing, Options::RADIAL_GRID_TYPES radialGridType,
                                              Options::SPHERICAL_GRID_TYPES sphericalGridType, unsigned int gridAccuracy,
                                              double gridWeightThreshold, const bool gridPointSorting) {
  /*
   * Store all the options into a grid factory
   */
  auto gridFactory = std::make_shared<GridFactory>(gridType, gridSmoothing, radialGridType, sphericalGridType,
                                                   gridAccuracy, gridWeightThreshold);
  /*
   * Construct and return the controller
   */
  return std::unique_ptr<AtomCenteredGridController>(new AtomCenteredGridController(geometry, gridFactory, gridPointSorting));
}

} /* namespace Serenity */
