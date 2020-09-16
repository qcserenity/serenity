/**
 * @file BasisExtension_test.cpp
 *
 * @date Feb 23, 2018
 * @author Moritz Bensberg
 * @copyright \n
 * This file is part of the program Serenity.\n\n
 * Serenity is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.\n\n
 * Serenity is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.\n\n
 * You should have received a copy of the GNU Lesser General
 * Public License along with Serenity.
 * If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Serenity Internal Headers */
#include "basis/BasisExtension.h"
#include "basis/BasisController.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class BasisExtensionTest
 * @brief Sets everything up for the tests of BasisExtension.h/.cpp .
 */
class BasisExtensionTest : public ::testing::Test {
 protected:
  BasisExtensionTest() {
  }

  virtual ~BasisExtensionTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(BasisExtensionTest, basisExtension) {
  auto sysA = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs);
  auto sysB = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs);
  BasisExtension basisExtender;
  basisExtender.extendAllSystems({sysA, sysB}, 1.0e-1);
  EXPECT_EQ((unsigned int)20, sysA->getBasisController()->getNBasisFunctions());
  EXPECT_EQ((unsigned int)24, sysB->getBasisController()->getNBasisFunctions());

  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
