/**
 * @file MolecularSurfaceFactory.h
 *
 * @author Moritz Bensberg
 * @date May 25, 2020
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

#ifndef GEOMETRY_MOLECULARSURFACEFACTORY_H_
#define GEOMETRY_MOLECULARSURFACEFACTORY_H_
/* Include Serenity Internal Headers */
#include "geometry/MolecularSurfaceController.h"
#include "misc/RememberingFactory.h"
#include "settings/PCMOptions.h"

namespace Serenity {

class PCMSettings;
class Atom;
class Geometry;

class MolecularSurfaceFactory
  : public RememberingFactory<MolecularSurfaceController, const std::shared_ptr<const Geometry>,
                              Options::PCM_CAVITY_TYPES, Options::PCM_ATOMIC_RADII_TYPES, unsigned int, double, double,
                              double, double, bool, unsigned int, unsigned int, double, double, bool, double> {
  MolecularSurfaceFactory() = default;

 public:
  virtual ~MolecularSurfaceFactory() = default;

  static std::shared_ptr<MolecularSurfaceController>
  produce(std::shared_ptr<const Geometry> geometry, Options::PCM_CAVITY_TYPES cavityType,
          Options::PCM_ATOMIC_RADII_TYPES radiiType, unsigned int calveLevel, double minDistance, double minRadius,
          double solvRadius, double overlapFactor, bool scaling, unsigned int sphericalAngularMomentum,
          unsigned int smallSphericalAngularMomentum, double alpha, double projectionCutOff, bool oneCavity,
          double connectivityFactor);

  static std::shared_ptr<MolecularSurfaceController> produce(std::shared_ptr<const Geometry> geometry,
                                                             const PCMSettings& pcmSettings);

 private:
  std::unique_ptr<MolecularSurfaceController>
  produceNew(std::shared_ptr<const Geometry> geometry, Options::PCM_CAVITY_TYPES cavityType,
             Options::PCM_ATOMIC_RADII_TYPES radiiType, unsigned int calveLevel, double minDistance, double minRadius,
             double solvRadius, double overlapFactor, bool scaling, unsigned int sphericalAngularMomentum,
             unsigned int smallSphericalAngularMomentum, double alpha, double projectionCutOff, bool oneCavity,
             double connectivityFactor) override final;

  double getAtomRadius(std::shared_ptr<Atom> atom, Options::PCM_ATOMIC_RADII_TYPES radiiType, bool scaling);
};

static std::unique_ptr<MolecularSurfaceFactory> _sufaceFactoryInstance;

} /* namespace Serenity */

#endif /* GEOMETRY_MOLECULARSURFACEFACTORY_H_ */
