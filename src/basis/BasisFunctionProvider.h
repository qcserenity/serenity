/**
 * @file   BasisFunctionProvider.h
 *
 * @date   24.03.2013
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
#ifndef BASISFUNCTIONPROVIDER_H_
#define BASISFUNCTIONPROVIDER_H_
/* Include Std and External Headers */
#include <string>

namespace Serenity {
class Atom;
/**
 * @class BasisFunctionProvider BasisFunctionProvider.h
 * @brief Provides basis functions using the Turbomole format
 *
 * as specified in the emsl basis set exchange library (https://bse.pnl.gov/bse/portal).
 */
class BasisFunctionProvider {
 private:
  /**
   * Private default constructor. Purely static class.
   */
  BasisFunctionProvider() = default;

 public:
  virtual ~BasisFunctionProvider() = default;
  /**
   * @brief Associates an atom with a set of @ref BasisFunction "BasisFunctions".
   *
   * @param atom      Will have (additional) basis functions assigned after call.
   * @param libraryPath A path to the folder in which the basis set files are stored.
   *                    The actual file name will be libraryPath+basisType.
   * @param basisType A label for the basis, e.g. 'def2-TZVP'.
   * @param isSpherical Whether the produced BasisFunctions shall be interpreted as spherical
   *                    (see BasisFunction class for details)
   * @param isPrimary Since an atom can have more than one set of basis functions this can be set
   *                  to mark whether the new set of basis functions shall be active (i.e. used per
   *                  default).
   */
  static void provideAtomWithBasisFunction(Atom& atom, const std::string libraryPath, const std::string basisType,
                                           bool isSpherical, bool isPrimary, int firstECP);

 private:
  static unsigned int resolveAngularMomentumChar(char type, std::string errorMessage);
};

} /* namespace Serenity */
#endif /* BASISFUNCTIONPROVIDER_H_ */
