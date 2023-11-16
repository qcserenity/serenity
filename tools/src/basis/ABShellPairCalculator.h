/**
 * @file ABShellPairCalculator.h
 *
 * @date Jun 20, 2018
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

#ifndef BASIS_ABSHELLPAIRCALCULATOR_H_
#define BASIS_ABSHELLPAIRCALCULATOR_H_
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class BasisController;
class ShellPairData;
/**
 * @class ABShellPairCalculator ABShellPairCalculator.h
 * @brief Calculates the ShellPairData for two arbitrary basis controllers.
 */
class ABShellPairCalculator {
 private:
  /*
   * Purely static.
   */
  ABShellPairCalculator();

 public:
  /**
   * @brief Calculates the ShellPairData for two basis controllers.
   * @param basisControllerA The basis controller A.
   * @param basisControllerB The basis controller B.
   * @return Returns the ShellPairData for the two basis controllers.
   */
  static std::shared_ptr<std::vector<ShellPairData>>
  calculateShellPairData_AB(std::shared_ptr<BasisController> basisControllerA, std::shared_ptr<BasisController> basisControllerB);
};

} /* namespace Serenity */

#endif /* BASIS_ABSHELLPAIRCALCULATOR_H_ */
