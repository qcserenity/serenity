/**
 * @file   GeometryAdderFactory.h
 *
 * @date   Dec 12, 2024
 * @author Anton Rikus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef GEOMETRYADDERFACTORY_H_
#define GEOMETRYADDERFACTORY_H_

/* Include Serenity Internal Headers */
#include "misc/RememberingFactory.h"
/* Include Std and External Headers */
#include <vector>

namespace Serenity {

/* Forward declarations */
class SystemController;
class Geometry;

/**
 * @class GeometryAdderFactory GeometryAdderFactory.h
 *
 * @brief The GeometryAdderFactory is a RememberingFactory for adding geometries of different systems. It was made since
 * multiple calls to the Kernel (where subsystem geometries were added) resulted in multiple constructions of the same
 * grid.
 *
 */
class GeometryAdderFactory
  : public RememberingFactory<Geometry, const std::vector<std::shared_ptr<SystemController>>, const std::vector<unsigned int>> {
 private:
  /**
   * Private default constructor - Singleton
   */
  GeometryAdderFactory() = default;

 public:
  virtual ~GeometryAdderFactory() = default;

  static std::shared_ptr<Geometry> produce(const std::vector<std::shared_ptr<SystemController>> activeSystems,
                                           const std::vector<unsigned int> subsystems);

 private:
  std::unique_ptr<Geometry> produceNew(const std::vector<std::shared_ptr<SystemController>> systems,
                                       const std::vector<unsigned int> systemIndices) override final;

  /*
   * Singleton: Instance of itself
   */
  static std::unique_ptr<GeometryAdderFactory> _instance;
};

} /* namespace Serenity */
#endif /* GEOMETRYADDERFACTORY_H_ */
