/**
 * @file   GridFactory.h
 *
 * @date   Mar 14, 2014
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
#ifndef GRIDFACTORY_H_
#define GRIDFACTORY_H_
/* Include Serenity Internal Headers */
#include "grid/AtomCenteredGrid.h"

namespace Serenity {
/* Forward declarations */
class AtomGridAccuracyDeterminator;
class AtomGridFactory;
class Geometry;
namespace Options {
enum class GRID_TYPES;
enum class RADIAL_GRID_TYPES;
enum class SPHERICAL_GRID_TYPES;
} // namespace Options

/**
 * @class GridFactory GridFactory.h
 *
 * @brief produces full integration grids, atom-centered and based on the atoms' integration grids
 */
class GridFactory {
 public:
  /**
   * @param flavour look into the Options.h for an explanation of the grid types.
   * @param smoothing in case of using the Becke partitioning
   * @param radialType for the construction of atom centered grids which are later combined
   * @param sphericalType for the construction of atom centered grids which are later combined
   * @param accuracyLevel simple and general accuracy parameter for the construction of the
   *                       atom-centered grids
   * @param weightThreshold grid points with a resulting weight smaller than this threshold will be
   *                        discarded
   */
  GridFactory(Options::GRID_TYPES flavour, unsigned int smoothing, const Options::RADIAL_GRID_TYPES radialType,
              const Options::SPHERICAL_GRID_TYPES sphericalType, unsigned int accuracyLevel, double weightThreshold);

  virtual ~GridFactory() = default;
  /**
   * @param geometry
   * @returns a new integration grid upon each call.
   */
  std::unique_ptr<AtomCenteredGrid> produce(std::shared_ptr<const Geometry> geometry);

 private:
  constexpr static double pRecursive(double nu, unsigned int k);
  constexpr static double beckeSmoothFunction(double nu, unsigned int k);

  Options::GRID_TYPES _flavour;
  unsigned int _smoothing;
  const Options::RADIAL_GRID_TYPES _radialType;
  const Options::SPHERICAL_GRID_TYPES _sphericalType;
  unsigned int _accuracyLevel;
  const double _weightThreshold;
};

} /* namespace Serenity */
#endif /* GRIDFACTORY_H_ */
