/**
 * @file   BasisController__TEST_SUPPLY.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   23. August 2015, 15:29
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
#ifndef BASISCONTROLLER__TEST_SUPPLY_H
#define BASISCONTROLLER__TEST_SUPPLY_H
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
/* Include Std and External Headers */
#include <map>
#include <memory>

namespace Serenity {
/**
 * These kinds of test bases are available:\n
 *
 * MINIMAL:\n
 *
 * A tiny basis with only two s-type functions.\n
 *
 * SMALL_MIXED:\n
 *
 * A small basis with an s-, p- and d-function. P-functions pose more complications than s-functions,
 * and d-functions pose more complications than p-functions. After that (i.e. f- and higher
 * functions) do not add conceptionally different complications. Thus with this basis already most
 * problems should be found.\n
 *
 * ASYMMETRIC_CH3:\n
 * \n
 * 6-31G* basis for an asymmetric CH3 radical\n
 *\n
 * 4\n
 *\n
 * C         -0.57388        1.74108        0.00000\n
 * H         -1.01199        0.28861        1.33493\n
 * H         -1.05807        0.62810        0.08722\n
 * H         -1.15616        2.73228        0.30975\n
 * \n
 *
 */
enum class TEST_BASIS_CONTROLLERS { MINIMAL, SMALL_MIXED, ASYMMETRIC_CH3 };
/**
 * @class BasisController__TEST_SUPPLY BasisController__TEST_SUPPLY.h
 * @brief Provides objects of type AtomCenteredBasisController without further dependencies.
 *
 * To be used as hard-coded test examples to work with in (unit) tests.
 */
class BasisController__TEST_SUPPLY {
  BasisController__TEST_SUPPLY() = delete;

 public:
  /**
   * @param kind specifies which test basis is requested
   * @returns a completely set up controller for the basis specified by kind
   */
  static std::shared_ptr<BasisController> getBasisController(TEST_BASIS_CONTROLLERS kind) {
    if (!_testBasisControllers[kind])
      prepare(kind);
    return std::dynamic_pointer_cast<BasisController>(_testBasisControllers[kind]);
  }
  /**
   * @param kind specifies which test basis is requested
   * @returns a completely set up controller for the basis specified by kind
   */
  static std::shared_ptr<AtomCenteredBasisController> getAtomCenteredBasisController(TEST_BASIS_CONTROLLERS kind) {
    if (!_testBasisControllers[kind])
      prepare(kind);
    return _testBasisControllers[kind];
  }

 private:
  static void prepare(TEST_BASIS_CONTROLLERS kind);
  static std::map<TEST_BASIS_CONTROLLERS, std::shared_ptr<AtomCenteredBasisController>> _testBasisControllers;
};

} /* namespace Serenity */
#endif /* BASISCONTROLLER__TEST_SUPPLY_H */
