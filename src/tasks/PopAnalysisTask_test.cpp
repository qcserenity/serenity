/**
 * @file PopAnalysisTask_test.cpp
 *
 * @date Aug 22, 2017
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "tasks/PopAnalysisTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <istream>

namespace Serenity {
/**
 * @class PopAnalysisTaskTest
 * @brief Sets everything up for the tests of PopAnalysisTask.h/.cpp .
 */
class PopAnalysisTaskTest : public ::testing::Test {
 protected:
  PopAnalysisTaskTest() {
  }

  virtual ~PopAnalysisTaskTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
* @brief Tests restricted output.
 */
TEST_F(PopAnalysisTaskTest, Mulliken_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  PopAnalysisTask<RESTRICTED> task(act);
  task.settings.mulliken = true;
  task.run();
}
/**
 * @test
 * @brief Tests unrestricted output.
 */
TEST_F(PopAnalysisTaskTest, Mulliken_U_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  act->setSCFMode(UNRESTRICTED);
  PopAnalysisTask<UNRESTRICTED> task(act);
  task.settings.mulliken = true;
  task.run();
}

} /*namespace Serenity*/
