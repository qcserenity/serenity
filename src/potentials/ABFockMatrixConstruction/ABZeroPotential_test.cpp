/**
 * @file ABZeroPotential_test.cpp
 *
 * @date Jun 26, 2018
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

/* Include Serenity Internal Headers */
#include "potentials/ABFockMatrixConstruction/ABZeroPotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABZeroPotentialTest : public ::testing::Test {
 protected:
  ABZeroPotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~ABZeroPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
};

/**
 * @test
 * @brief Tests the restricted construction of the f_AB matrix.
 */
TEST_F(ABZeroPotentialTest, H2_dimer_rFockMatrixABC) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();

  ABZeroPotential<SPIN> zeroPot(basisA, basisB);

  SPMatrix<SPIN> f_AB = zeroPot.getMatrix();

  EXPECT_NEAR(f_AB(0, 0), 0.0, 1e-10);
  EXPECT_NEAR(f_AB(1, 0), 0.0, 1e-10);
  EXPECT_NEAR(f_AB(2, 0), 0.0, 1e-10);
  EXPECT_NEAR(f_AB(0, 1), 0.0, 1e-10);
  EXPECT_NEAR(f_AB(0, 2), 0.0, 1e-10);
  EXPECT_NEAR(f_AB(0, 3), 0.0, 1e-10);
}

} /* namespace Serenity */
