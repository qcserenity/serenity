/**
 * @file   CombinedShellPair.h
 *
 * @date   Jul 30, 2020
 * @author Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef SRC_BASIS_COMBINEDSHELLPAIR_H_
#define SRC_BASIS_COMBINEDSHELLPAIR_H_

/* Include Serenity Internal Headers */
#include "basis/Shell.h"
#include "notification/NotifyingClass.h"

namespace Serenity {

/**
 * @class CombinedShellPair CombinedShellPair.h
 *
 * A class that holds the product of two shells as a regular shell.
 * Most importantly, this class assures that the correct normalization
 * of the base shells is used so that \f$ ( i j | k l ) = ( [ij] | [kl] ) \f$.
 *
 */
class CombinedShellPair : public Shell {
 public:
  /**
   * @brief Constructor
   *
   * @param shellA The first shell in the shell product.
   * @param shellB The second shell in the shell product.
   * @param spherical A Flag if spherical or cartesian basis functions are used.
   */
  CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB, bool spherical);

  /**
   * @brief Constructor that allows to manually set the angular momentum
   *                of the shell product. This also adjusts the normalization accordingly.
   *
   * @param shellA The first shell in the shell product.
   * @param shellB The second shell in the shell product.
   * @param angularMomentum The angular momentum of the shell product.
   * @param spherical A Flag if spherical or cartesian basis functions are used.
   */
  CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                    unsigned int angularMomentum, bool spherical);

  /**
   * @brief Constructor that allows to manually set the angular momentum, the exponents
   *                and contractions  of the shell product. This also adjusts the normalization accordingly.
   *
   * @param shellA The first shell in the shell product.
   * @param shellB The second shell in the shell product.
   * @param expo The exponents used for the combined shell.
   * @param contr The contractions for the combined shells.
   * @param angularMomentum The angular momentum of the shell product.
   * @param spherical A Flag if spherical or cartesian basis functions are used.
   */
  CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB, std::vector<double> expo,
                    std::vector<double> contr, unsigned int angularMom, bool spherical);
  /**
   * @brief Constructor that allows to manually set the angular momentum, the exponents
   *                and contractions  of the shell product. This also adjusts the normalization accordingly.
   *
   * @param shellA The first shell in the shell product.
   * @param shellB The second shell in the shell product.
   * @param expo The exponents used for the combined shell.
   * @param contr The contractions for the combined shells.
   * @param angularMomentum The angular momentum of the shell product.
   * @param spherical A Flag if spherical or cartesian basis functions are used.
   * @param coords The coordinates of the base atom.
   * @param element The element string of the base atom.
   */
  CombinedShellPair(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB, std::vector<double> expo,
                    std::vector<double> contr, unsigned int angularMomentum, bool spherical,
                    std::array<double, 3> coords, std::string element);

  /**
   * @brief A getter for the first base shell of the product
   *
   * @return The first base shell.
   */
  std::shared_ptr<const Shell> getBaseShellA() {
    return _shellA;
  }
  /**
   * @brief A getter for the second base shell of the product
   *
   * @return The second base shell.
   */
  std::shared_ptr<const Shell> getBaseShellB() {
    return _shellB;
  }

 private:
  /**
   * @brief Generates the set of unique exponents for the shell product
   *
   * @param shellA The first shell to get exponents from.
   * @param shellB The second shell to get exponents from.
   * @return
   */
  std::vector<double> generateExponents(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB);

  /**
   * @brief Generates the set of unique correctly normalized contraction
   *                coefficients for the shell product.
   *
   * @param shellA The first shell to get contractions from.
   * @param shellB The second shell to get contractions from.
   * @param angularMomentum The angular momentum assigned to the shell contractions are generated for.
   * @return
   */
  std::vector<double> generateContractions(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB,
                                           unsigned int angularMomentum);
  /**
   * @brief Generates the set of unique contraction
   *                coefficients for the shell product without additional renormalization.
   *
   * @param shellA The first shell to get contractions from.
   * @param shellB The second shell to get contractions from.
   * @return
   */
  std::vector<double> generateContractionsNoNorm(std::shared_ptr<const Shell> shellA, std::shared_ptr<const Shell> shellB);

  /**
   * @brief Reverses the normalization that is always done by libint. This is important as the
   *                combined shell should represent the normalized base shells and not the normalized
   *                combined shell.
   * @param expo The exponents of the shell.
   * @param contr The contractions of the shell.
   * @param angularMomentum The angular momentum of the shell.
   * @return The coefficients scalled by the inverse of their normalization factors.
   */
  std::vector<double> reverseNormalization(std::vector<double> expo, std::vector<double> contr, unsigned int angularMomentum);
  /**
   * @brief Function used to assert that both shells are defined on the same center
   *                during construction.
   *
   * @param a The coordinates for the first atom.
   * @param b The coordinates for the second atom.
   * @return The coordinates of the first system if they are equal to those of the second.
   */
  std::array<double, 3> checkCoords(std::array<double, 3> a, std::array<double, 3> b);
  /**
   * @brief Function used to assert that both shells are defined for the same element
   *                during construction.
   *
   * @param a The element string for the first atom.
   * @param b The element string for the second atom.
   * @return The element string of the first system if it is equal to those of the second.
   */
  std::string checkElement(std::string a, std::string b);
  /// The first base shell.
  std::shared_ptr<const Shell> _shellA;
  /// The second base shell.
  std::shared_ptr<const Shell> _shellB;
};

} // namespace Serenity
#endif /* SRC_BASIS_COMBINEDSHELLPAIR_H_ */
