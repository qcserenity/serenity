/**
 * @file AtomCenteredBasisControllerFactory.h
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
#ifndef ATOMCENTEREDBASISCONTROLLERFACTORY_H
#define ATOMCENTEREDBASISCONTROLLERFACTORY_H
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "misc/RememberingFactory.h"
/* Include Std and External Headers */
#include <memory>
#include <string>

namespace Serenity {
/* Forward Declarations */
class Geometry;
/**
 * @class AtomCenteredBasisControllerFactory AtomCenteredBasisControllerFactory.h
 * @brief produces instances of AtomCenteredBasisController
 */
class AtomCenteredBasisControllerFactory
  : public RememberingFactory<AtomCenteredBasisController, const std::shared_ptr<Geometry>, const std::string, bool,
                              bool, const std::string, int> {
 private:
  AtomCenteredBasisControllerFactory() = default;

 public:
  virtual ~AtomCenteredBasisControllerFactory() = default;
  /**
   * @param geometry
   * @param basisLibrary          Path to the basis set library (folder)
   * @param basisLabel            The name of the basis e.g. "6-31Gs"
   * @param makeSphericalBasis    Whether or not to transform the basis in the spherical space.
   * @param makePrimary           This switch will set the basis as the primary one in each atom.
   * @param firstECP              Nuclear charge threshold at which ECPs are used if available.
   * @returns a new or already existing AtomCenteredBasisController
   */
  static std::shared_ptr<AtomCenteredBasisController>
  produce(const std::shared_ptr<Geometry> geometry, const std::string basisLibrary, bool makeSphericalBasis,
          bool makePrimary, int firstECP, const std::string basisLabel = "");

 private:
  std::unique_ptr<AtomCenteredBasisController>
  produceNew(const std::shared_ptr<Geometry> geometry, const std::string basisLibrary, bool makeSphericalBasis,
             bool makePrimary, const std::string basisLabel, int firstECP) override final;
  /*
   * Singleton: Instance of itself
   */
  static std::unique_ptr<AtomCenteredBasisControllerFactory> _instance;
};

} /* namespace Serenity */
#endif /* ATOMCENTEREDBASISCONTROLLERFACTORY_H */
