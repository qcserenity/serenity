/**
 * @file   AtomGridFactory.h
 *
 * @date   Mar 14, 2014
 * @author Dennis Barton
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
#ifndef ATOMGRIDFACTORY_H
#define ATOMGRIDFACTORY_H
/* Include Serenity Internal Headers */
#include "misc/RememberingFactory.h"
#include "parameters/AtomicParameters.h"
/* Include Std and External Headers */
#include <array>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
namespace Options {
enum class RADIAL_GRID_TYPES;
enum class SPHERICAL_GRID_TYPES;
} // namespace Options
class AtomGrid;
class AtomType;
/**
 * @class AtomGridFactory AtomGridFactory.h
 *
 * @brief Construction of atom specific grid points
 *
 * Construct grid points based of atom type from angular
 * points (mostly Lebedev others possible) and radial
 * points from specific Gauss-Chebychev (Becke, Handy, Ahlrichs etc.)
 */
class AtomGridFactory
  : public RememberingFactory<const AtomGrid, const Options::RADIAL_GRID_TYPES, const Options::SPHERICAL_GRID_TYPES,
                              const std::shared_ptr<const AtomType>, const unsigned int> {
 private:
  /**
   * Private default constructor - Singleton
   */
  AtomGridFactory() = default;

 public:
  virtual ~AtomGridFactory() = default;
  /**
   * @brief Here are the actual atom centric points constructed
   *
   * construct grid based on the following parameters
   * @param radType radial grid type
   * @param sphType spherical grid type
   * @param atomType Type of Atom to construct a grid for
   * @param accuracy The accuracy [1,7]
   * @returns an atom specific grid
   */
  static std::shared_ptr<const AtomGrid> produce(Options::RADIAL_GRID_TYPES radType, Options::SPHERICAL_GRID_TYPES sphType,
                                                 std::shared_ptr<const AtomType> atomType, unsigned int accuracy);
  /**
   * @param radType radial grid type
   * @param sphType spherical grid type
   * @param atomType Type of Atom to construct a grid for
   * @param accuracy The accuracy [1,7]
   * @returns an atom specific grid
   */
  std::unique_ptr<const AtomGrid>
  produceNew(const Options::RADIAL_GRID_TYPES radialType, const Options::SPHERICAL_GRID_TYPES sphericalType,
             const std::shared_ptr<const AtomType> atomType, const unsigned int accuracy) override final;

  /*
   * The following functions should actually be private, but are directly used during testing and thus
   * set to public.
   */
  /**
   * @brief Here the radial points are calculated
   *
   * Depending of the value of _radType
   * we'll construct the formula for radial integration
   * corresponding to the give type
   *
   * @param alpha atom dependent scaling factor for abscissa
   * @param nRadial number of radial grid points
   * @param radPoints vector of radial grid points
   * @param radWeights integration weights for corresponding grid point
   * @param radType see RADIAL_GRID_TYPES
   */
  static void _radialGrid(double alpha, unsigned int nRadial, std::vector<double>& radPoints,
                          std::vector<double>& radWeights, const Options::RADIAL_GRID_TYPES& radType);
  /**
   * @brief Here lebedev grid points are calculated
   *
   * Depending on the i value (which should be in the range of
   * 0 to approx 30 to 40) directly connected to the angular order
   * of the Lebedev grid i.e.
   * i=0 should give 6 angular points, i=1 should give 14,
   * and so on...\n
   *
   * Look into the sphere_lebedev_rule.cpp for a reference.
   *
   * @param[in] i value corresponding to the angular order (i.e. determines granularity)
   * @param[out] sphN determines the size of the following vectors, i.e. the number of points
   * @param[out] sphX x-value of spherical grid point
   * @param[out] sphY y-value of    "        "    "
   * @param[out] sphZ z-value of    "        "    "
   * @param[out] sphW weight of corresponding spherical grid point
   */
  static void _lebedevSphericalGrid(unsigned int i, unsigned int& sphN, std::vector<double>& sphX,
                                    std::vector<double>& sphY, std::vector<double>& sphZ, std::vector<double>& sphW);

 private:
  /*
   * Values for the atomic size modifications as suggested by
   * Treutler, Ahlrichs, J. Chem. Phys. 102 (1995), 346.
   */
  static const std::array<double, N_ELEMENTS_IN_PERIODIC_TABLE + 1> _ahlrichsAlphaValues;
  /*
   * Calculated atomic radii as proposed by
   * E. Clementi, D. L. Raimondi, W. P. Reinhardt, J. Chem. Phys. 47 (1967), 1300
   */
  static const std::array<double, N_ELEMENTS_IN_PERIODIC_TABLE + 1> _clementisRadii;
  static std::unique_ptr<AtomGridFactory> _instance;
};

} /* namespace Serenity */

#endif /* ATOMGRIDFACTORY_H */
