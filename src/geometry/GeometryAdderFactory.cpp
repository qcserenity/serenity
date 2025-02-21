/**
 * @file GeometryAdderFactory.cpp
 *
 * @date Dec 12, 2024
 * @author Anton Rikus
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
#include "geometry/GeometryAdderFactory.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "system/SystemController.h"

namespace Serenity {

/*
 * Initialize static pointer for the singleton instance.
 */
std::unique_ptr<GeometryAdderFactory> GeometryAdderFactory::_instance;

std::shared_ptr<Geometry> GeometryAdderFactory::produce(const std::vector<std::shared_ptr<SystemController>> systems,
                                                        const std::vector<unsigned int> systemIndices) {
  if (!_instance)
    _instance.reset(new GeometryAdderFactory());
  return _instance->getOrProduce(systems, systemIndices);
}

std::unique_ptr<Geometry> GeometryAdderFactory::produceNew(const std::vector<std::shared_ptr<SystemController>> systems,
                                                           const std::vector<unsigned int> systemIndices) {
  std::unique_ptr<Geometry> geometry = std::make_unique<Geometry>();

  // Build geometry from subsystem geometries.
  if (systemIndices.size() != 0) {
    for (unsigned I = 0; I < systemIndices.size(); ++I) {
      if (systemIndices[I] > systems.size()) {
        throw SerenityError("Your subsystem index is larger than the number of subsystems!");
      }
      (*geometry) += (*systems[systemIndices[I] - 1]->getGeometry());
    }
  }
  // Build supersystem grid.
  else {
    if (systems.size() == 1) {
      *geometry += *systems[0]->getGeometry();
    }
    else {
      for (unsigned I = 0; I < systems.size(); ++I) {
        (*geometry) += (*systems[I]->getGeometry());
      }
    }
  }
  geometry->deleteIdenticalAtoms();
  return geometry;
}

} /* namespace Serenity */
