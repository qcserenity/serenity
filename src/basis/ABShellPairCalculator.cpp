/**
 * @file ABShellPairCalculator.cpp
 *
 * @date   Dec 16, 2020
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
/* Include Class Header*/
#include "basis/ABShellPairCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h" //Loop shells.
#include "basis/BasisController.h"
#include "integrals/wrappers/Libint.h"

namespace Serenity {

std::shared_ptr<std::vector<ShellPairData>>
ABShellPairCalculator::calculateShellPairData_AB(std::shared_ptr<BasisController> basisControllerA,
                                                 std::shared_ptr<BasisController> basisControllerB) {
  // intialize libint
  auto& libint = Libint::getInstance();
  libint.initialize_plain(LIBINT_OPERATOR::coulomb, 4, std::numeric_limits<double>::epsilon() / 1e4,
                          std::max(basisControllerA->getMaxNumberOfPrimitives(), basisControllerB->getMaxNumberOfPrimitives()));
  auto shellPairList = std::make_shared<std::vector<ShellPairData>>();
  const auto& basisA = basisControllerA->getBasis();
  const auto& basisB = basisControllerB->getBasis();
  // loops over shells
  Eigen::MatrixXd integrals;
  const unsigned int nShellsA = basisA.size();
  const unsigned int nShellsB = basisB.size();
  for (unsigned int a = 0; a < nShellsA; ++a) {
    const auto& shellA = *basisA[a];
    for (unsigned int b = 0; b < nShellsB; ++b) {
      const auto& shellB = *basisB[b];
      // calculate integrals
      if (libint.compute(LIBINT_OPERATOR::coulomb, 0, shellA, shellB, shellA, shellB, integrals)) {
        (*shellPairList).push_back(ShellPairData(a, b, sqrt(integrals.maxCoeff()), false));
      } /* if (prescreen) */
    }   /* b/shellB */
  }     /* a/shellA */
  // finalize libint
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 4);
  // sort the list
  std::sort((*shellPairList).begin(), (*shellPairList).end());
  std::reverse((*shellPairList).begin(), (*shellPairList).end());
  return shellPairList;
};

} // namespace Serenity
