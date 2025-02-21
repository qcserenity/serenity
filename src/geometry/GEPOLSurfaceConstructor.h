/**
 * @file GEPOLSurfaceConstructor.h
 *
 * @date Sep 27, 2016
 * @author David Schnieders
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

#ifndef GEOMETRY_GEPOLSURFACECONSTRUCTOR_H_
#define GEOMETRY_GEPOLSURFACECONSTRUCTOR_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
class Sphere;
class Triangle;
class Atom;
class Point;
class GridController;
class MolecularSurface;
struct PCMSettings;
namespace Options {
enum class PCM_ATOMIC_RADII_TYPES;
enum class PCM_CAVITY_TYPES;
} // namespace Options
/**
 * @class GEPOLSurfaceConstructor GEPOLSurfaceConstructor.h
 * @brief Constructs the GEPOL surface.
 */
class GEPOLSurfaceConstructor {
 public:
  /**
   * @brief Constructor.
   * @param spheres           The initial spheres.
   * @param isSAS             Construct SAS surface.
   * @param patchLevel        Level of increased resolution at sphere intersections.
   * @param minDistance       Minimal distance between triangle centeres.
   * @param minRadius         Minimal radius of added spheres.
   * @param probeRadius       Radius of the solvent probe.
   * @param overlapFactor     Maximum overlap between spheres.
   */
  GEPOLSurfaceConstructor(std::vector<Sphere> spheres, bool isSAS, unsigned int patchLevel, double minDistance,
                          double minRadius, double probeRadius, double overlapFactor);
  /**
   * @brief Default destructor.
   */
  virtual ~GEPOLSurfaceConstructor();
  /**
   * @brief Getter for the final molecular surface
   * @return The molecular surface.
   */
  std::unique_ptr<MolecularSurface> getMolecularSurface();

 private:
  // Construct the surface from spheres.
  void buildSurface(const std::vector<Sphere>& spheres);
  // Initialize the surface.
  void initializeSurface();
  // Construct all spheres.
  void addSpheres(bool coarse, std::vector<Sphere>& spheres);
  // Search for overlapping triangles and covered triangle centers.
  void checkTriangleCenter(bool& withinSurface, bool& intersected, const Triangle& triangle, const double& effSolvRad,
                           const std::vector<Sphere>& spheres);
  // Increase the resolution at an intersection.
  std::vector<Triangle> patchIntersection(const Triangle& triangle, double effSolvRad, int patchLevel,
                                          const std::vector<Sphere>& spheres);
  // If true, the surface is available.
  bool _initialized;
  // Flag for solvent accessible surface (SAS).
  bool _isSAS;
  // Reselution increase at intersection.
  const unsigned int _calveLevel;
  // Minimal distance between triangle centers.
  double _minDistance;
  // Minimal radius of added spheres.
  double _minRad;
  // Radius of the solvent probe.
  double _solvRad;
  // Maximum overlap between spheres.
  double _overlapFactor;
  // The spheres.
  std::vector<Sphere> _spheres;
  // The triangles.
  std::vector<Triangle> _molecularSurface;
  // The molecular surface..
  std::unique_ptr<MolecularSurface> _surface;
};

} /* namespace Serenity */

#endif /* GEOMETRY_GEPOLSURFACECONSTRUCTOR_H_ */
