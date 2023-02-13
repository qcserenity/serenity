/**
 * @file CartesianToSphericalTransformer.h
 *
 * @date May 24, 2018
 * @author Moritz Bensberg
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
#ifndef BASIS_CARTESIANTOSPHERICALTRANSFORMER_H_
#define BASIS_CARTESIANTOSPHERICALTRANSFORMER_H_
/* Include Serenity Internal Headers */
#include "math/IntegerMaths.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class CartesianToSphericalTransformer CartesianToSphericalTransformer.h
 * @brief This purely static class calculates transformation matrices for the transformation between
 *        spherical and Cartesian solid harmonics. The equations for the transformation can be found
 *        in "Molecular Electronic-Structure Theory", Helgaker, Jørgensen, Olsen, p. 215
 *        eq.6.4.47 to eq. 6.4.50.
 */
class CartesianToSphericalTransformer {
 public:
  /// @brief Default destructor.
  virtual ~CartesianToSphericalTransformer();
  /**
   * @brief Calculates the transformation matrix for l or gets it from memory.
   * @param l The angular momentum of the shell.
   * @return The transformation matrix.
   */
  static Eigen::MatrixXd& getTransformationMatrix(unsigned int l);

 private:
  // Purely static, never instantiated.
  CartesianToSphericalTransformer();

  /* Worker functions */
  // The normalization factor N^S_lm
  // Implements eq. 6.4.49 p. 215 Helgaker, Jørgensen, Olsen
  static double norm(unsigned int l, unsigned int abs_m) {
    return sqrt(2 * factorial(l + abs_m) * factorial(l - abs_m) / ((abs_m == 0) ? 2 : 1)) / (intPow(2, abs_m) * factorial(l));
  }
  // The transformation coefficient C^(lm)_(tuv).
  // Implements eq.  6.4.48 p. 215 Helgaker, Jørgensen, Olsen
  static double coef(unsigned int l, unsigned int abs_m, unsigned int t, unsigned int u, unsigned int vtimes2);
  // Maps between tuv and the index of the associated Cartesian harmonic. See exponents of x,y,z
  // in eq. 6.4.47 p. 215 Helgaker, Jørgensen, Olsen
  static unsigned int mapToCartHarmonics(unsigned int l, unsigned int abs_m, unsigned int ttimes2, unsigned int utimes2,
                                         unsigned int vtimes2);
  // The transformation matrices are calculated once and then stored in this map for the case
  // that they are used again. The matrices are fairly small so this should not result in memory issues.
  // Because this class is purely static the map has to be pre-defined. Therefore, only angular moments
  // l<=20 are covered in the map. This should be sufficient for pretty much any basis set.
  static std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> _transformationMatrices;
};

} /* namespace Serenity */

#endif /* BASIS_CARTESIANTOSPHERICALTRANSFORMER_H_ */
