/**
 * @file   GridControllerFactory.h
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
#ifndef GRIDCONTROLLERFACTORY_H_
#define GRIDCONTROLLERFACTORY_H_
/* Include Serenity Internal Headers */
#include "grid/AtomCenteredGridController.h"
#include "misc/RememberingFactory.h"
#include "settings/GridOptions.h"

namespace Serenity {
/* Forward declarations */
class Geometry;
class GridFactory;
struct Settings;
struct GRID;

/**
 * TODO rename to AtomCenteredGridControllerFactory
 * @class GridControllerFactory GridControllerFactory.h
 * @brief Produces instances of GridController
 */
class GridControllerFactory
  : public RememberingFactory<AtomCenteredGridController, const std::shared_ptr<const Geometry>, Options::GRID_TYPES, unsigned int,
                              Options::RADIAL_GRID_TYPES, Options::SPHERICAL_GRID_TYPES, unsigned int, double, const bool> {
 private:
  /**
   * Private default constructor, singleton.
   */
  GridControllerFactory() = default;

 public:
  virtual ~GridControllerFactory() = default;
  /**
   * @param geometry
   * @returns a new or already existing GridController
   */
  static std::shared_ptr<AtomCenteredGridController>
  produce(std::shared_ptr<const Geometry> geometry, Options::GRID_TYPES gridType, unsigned int gridSmoothing,
          Options::RADIAL_GRID_TYPES radialGridType, Options::SPHERICAL_GRID_TYPES sphericalGridType,
          unsigned int gridAccuracy, double gridWeightThreshold, const bool gridPointSorting);
  /**
   * @param geometry
   * @param settings from here the parameters for the grids the controller will work with are taken
   * @param gridPurpose to choose between regular and small grid
   * @returns a new or already existing GridController
   */
  static std::shared_ptr<AtomCenteredGridController> produce(std::shared_ptr<const Geometry> geometry, const GRID& gridSettings,
                                                             Options::GRID_PURPOSES purpose = Options::GRID_PURPOSES::DEFAULT);

 private:
  std::unique_ptr<AtomCenteredGridController>
  produceNew(const std::shared_ptr<const Geometry> geometry, Options::GRID_TYPES gridType, unsigned int gridSmoothing,
             Options::RADIAL_GRID_TYPES radialGridType, Options::SPHERICAL_GRID_TYPES sphericalGridType,
             unsigned int gridAccuracy, double gridWeightThreshold, const bool gridPointSorting) override final;
};

static std::unique_ptr<GridControllerFactory> _instance;

} /* namespace Serenity */

#endif /* GRIDCONTROLLERFACTORY_H_ */
