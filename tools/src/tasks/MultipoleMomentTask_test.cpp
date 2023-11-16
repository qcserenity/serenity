/**
 * @file MultipoleMomentTask_test.cpp
 *
 * @date Aug 22, 2017
 * @author Jan Unsleber
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
#include "tasks/MultipoleMomentTask.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <istream>

namespace Serenity {
/**
 * @class MultipoleMomentTaskTest
 * @brief Sets everything up for the tests of MultipoleMomentTask.h/.cpp .
 */
class MultipoleMomentTaskTest : public ::testing::Test {
 protected:
  MultipoleMomentTaskTest() {
  }

  virtual ~MultipoleMomentTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests restricted output.
 */
TEST_F(MultipoleMomentTaskTest, Quad_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  MultipoleMomentTask task({act});
  task.settings.highestOrder = 2;
  task.run();
}
/**
 * @test
 * @brief Tests unrestricted output.
 */
TEST_F(MultipoleMomentTaskTest, Quad_U_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  act->setSCFMode(UNRESTRICTED);
  MultipoleMomentTask task({act});
  task.settings.highestOrder = 2;
  task.run();
}

} /*namespace Serenity*/
