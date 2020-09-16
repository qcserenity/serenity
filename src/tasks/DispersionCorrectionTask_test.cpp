/**
 * @file DispersionCorrectionTask_test.cpp
 *
 * @date    Nov 9, 2017
 * @author: Jan Unsleber
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
#include "tasks/DispersionCorrectionTask.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class DispersionCorrectionTaskTest
 * @brief Sets everything up for the tests of DispersionCorrectionTask.h/.cpp .
 */
class DispersionCorrectionTaskTest : public ::testing::Test {
 protected:
  DispersionCorrectionTaskTest() {
  }

  virtual ~DispersionCorrectionTaskTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests DispersionCorrectionTask.h/.cpp:
 *   Test if output is running without crashes for all corrections.
 */
TEST_F(DispersionCorrectionTaskTest, All) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto task = DispersionCorrectionTask(systemController);
  task.settings.gradient = true;
  task.settings.hessian = true;
  task.run();
}

} /* namespace Serenity */
