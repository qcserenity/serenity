/**
 * @file MP2_test.cpp
 *
 *  @date      Nov 24, 2015
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
#include "postHF/MPn/MP2.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class MP2Test
 * @brief Sets everything up for the tests of MP2.h/.cpp .
 */
class MP2Test : public ::testing::Test {
 protected:
  MP2Test() : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS)) {
    systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  }

  virtual ~MP2Test() = default;

  /// system
  std::shared_ptr<SystemController> systemController;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests MP2.h/.cpp: Restricted energy test.
 */
TEST_F(MP2Test, MP2Restricted) {
  MP2EnergyCorrector<Options::SCF_MODES::RESTRICTED> mp2EnergyCorrector(systemController);
  auto MP2EnergyCorrection = mp2EnergyCorrector.calculateElectronicEnergy();
  EXPECT_NEAR(MP2EnergyCorrection, -0.01407063069504663, 1E-8);
}

} // namespace Serenity
