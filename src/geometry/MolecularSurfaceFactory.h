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
#include "geometry/MolecularSurface.h"
#include "settings/PCMOptions.h"

namespace Serenity {

struct PCMSettings;
class Atom;
class Geometry;
/**
 * @class
 * @brief A factory to construct a molecular surface according to the given settings.
 */
class MolecularSurfaceFactory {
  MolecularSurfaceFactory() = default;

 public:
  virtual ~MolecularSurfaceFactory() = default;
  /**
   * @brief Produce a new molecular surface.
   * @param geometry                         The molecular geometry.
   * @param cavityType                       The type of cavity.
   * @param radiiType                        The radii set.
   * @param calveLevel                       Level of resolution increase at sphere intersections (GEPOL).
   * @param minDistance                      Minimum tolerated distance between surface points.
   * @param minRadius                        Minimum radius of any added sphere (GEPOL)
   * @param solvRadius                       Radius of the solvent probe.
   * @param overlapFactor                    Overlap ratio allowed for added spheres.
   * @param scaling                          If true, atom radii are scaled by 1.2.
   * @param sphericalAngularMomentum         Angular momentum for the Lebedev-spherical grids (DELLEY).
   * @param smallSphericalAngularMomentum    Angular momentum for the Lebedev-spherical grids for H-atoms (DELLEY).
   * @param alpha                            Sharpness paramenter for the DELLEY-surface.
   * @param projectionCutOff                 Projection cut-off for the DELLEY-surface.
   * @param oneCavity                        Remove cavity points that are not connected to the point with the most
   *                                         extreme x-coordinate (DELLEY).
   * @param connectivityFactor               Connection between grid points is established as a range cut-off. The
   * threshold is given by connectivityFactor times solvRadius.
   */
  static std::unique_ptr<MolecularSurface>
  produce(std::shared_ptr<const Geometry> geometry, Options::PCM_CAVITY_TYPES cavityType,
          Options::PCM_ATOMIC_RADII_TYPES radiiType, unsigned int calveLevel, double minDistance, double minRadius,
          double solvRadius, double overlapFactor, bool scaling, unsigned int sphericalAngularMomentum,
          unsigned int smallSphericalAngularMomentum, double alpha, double projectionCutOff, bool oneCavity,
          double connectivityFactor);
  /**
   * @brief Produce a new molecular surface.
   * @param The geometry.
   * @param The PCM settings.
   */
  static std::unique_ptr<MolecularSurface> produce(std::shared_ptr<const Geometry> geometry, const PCMSettings& pcmSettings);

 private:
  // Helper function to select the right atom radius.
  static double getAtomRadius(std::shared_ptr<Atom> atom, Options::PCM_ATOMIC_RADII_TYPES radiiType, bool scaling);
};

} /* namespace Serenity */

#endif /* GEOMETRY_MOLECULARSURFACEFACTORY_H_ */
