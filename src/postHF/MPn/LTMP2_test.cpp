/**
 * @file LTMP2_test.cpp
 *
 * @date Dec 30, 2022
 * @author Niklas Niemeyer
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

#ifdef SERENITY_USE_LAPLACE_MINIMAX

/* Include Serenity Internal Headers */
#include "postHF/MPn/LTMP2.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LTMP2Test
 * @brief Sets everything up for the tests of LTMP2.h/.cpp .
 */
class LTMP2Test : public ::testing::Test {
 protected:
  LTMP2Test()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_HF)) {
    systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  }

  virtual ~LTMP2Test() = default;

  /// system
  std::shared_ptr<SystemController> systemController;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests LTMP2.h/.cpp: Restricted energy test.
 */
TEST_F(LTMP2Test, MP2Restricted) {
  LTMP2<Options::SCF_MODES::RESTRICTED> LTMP2(systemController);
  auto MP2EnergyCorrection = LTMP2.calculateCorrection();
  EXPECT_NEAR(MP2EnergyCorrection, -0.036065315338922259, 1E-8);
}

/**
 * @test
 * @brief Tests LTMP2.h/.cpp: Unrestricted energy test.
 */
TEST_F(LTMP2Test, MP2Unrestricted) {
  LTMP2<Options::SCF_MODES::UNRESTRICTED> LTMP2(systemController);
  auto MP2EnergyCorrection = LTMP2.calculateCorrection();
  EXPECT_NEAR(MP2EnergyCorrection, -0.036065315338922259, 1E-8);
}

} // namespace Serenity

#endif /* SERENITY_USE_LAPLACE_MINIMAX */