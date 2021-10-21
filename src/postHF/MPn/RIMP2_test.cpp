/**
 * @file RIMP2_test.cpp
 *
 *  @date      Aug 14, 2017
 *  @author    Jan Unsleber
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
#include "postHF/MPn/RIMP2.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class RIMP2Test
 * @brief Sets everything up for the tests of RIMP2.h/.cpp .
 */
class RIMP2Test : public ::testing::Test {
 protected:
  RIMP2Test()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS)) {
    systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  }

  virtual ~RIMP2Test() = default;

  /// system
  std::shared_ptr<SystemController> systemController;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests RIMP2.h/.cpp: Restricted energy test.
 */
TEST_F(RIMP2Test, MP2Restricted) {
  RIMP2<Options::SCF_MODES::RESTRICTED> rimp2(systemController);
  auto MP2EnergyCorrection = rimp2.calculateCorrection();
  EXPECT_NEAR(MP2EnergyCorrection, -0.014069347723441622, 1E-8);
}

/**
 * @test
 * @brief Tests RIMP2.h/.cpp: Unrestricted energy test.
 */
TEST_F(RIMP2Test, MP2Unrestricted) {
  RIMP2<Options::SCF_MODES::UNRESTRICTED> rimp2(systemController);
  auto MP2EnergyCorrection = rimp2.calculateCorrection();
  EXPECT_NEAR(MP2EnergyCorrection, -0.014069347723441622, 1E-8);
}

} // namespace Serenity
