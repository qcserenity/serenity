/**
 * @file AtomCenteredBasisControllerFactory.cpp
 *
 * @date Jul 31, 2015
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
#include "basis/AtomCenteredBasisControllerFactory.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

/*
 * Initialize static pointer for the singleton instance.
 */
std::unique_ptr<AtomCenteredBasisControllerFactory> AtomCenteredBasisControllerFactory::_instance;

std::shared_ptr<AtomCenteredBasisController>
AtomCenteredBasisControllerFactory::produce(std::shared_ptr<Geometry> geometry, const std::string basisLibrary,
                                            bool makeSphericalBasis, bool makePrimary, int firstECP,
                                            const std::string basisLabel) {
  assert(geometry);
  std::string copy = basisLabel;
  for (auto& c : copy)
    c = std::toupper(c);
  if (!_instance)
    _instance.reset(new AtomCenteredBasisControllerFactory());
  return _instance->getOrProduce(geometry, basisLibrary, makeSphericalBasis, makePrimary, copy, firstECP);
}

std::unique_ptr<AtomCenteredBasisController>
AtomCenteredBasisControllerFactory::produceNew(const std::shared_ptr<Geometry> geometry, const std::string basisLibrary,
                                               bool makeSphericalBasis, bool makePrimary, const std::string basisLabel,
                                               int firstECP) {
  return std::unique_ptr<AtomCenteredBasisController>(
      new AtomCenteredBasisController(geometry, basisLibrary, makeSphericalBasis, makePrimary, basisLabel, firstECP));
}

} /* namespace Serenity */
