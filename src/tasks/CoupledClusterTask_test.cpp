/**
 * @file CoupledClusterTask_test.cpp
 *
 * @date Nov 17, 2017
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
#include "tasks/CoupledClusterTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <istream>

namespace Serenity {
/**
 * @class PopAnalysisTaskTest
 * @brief Sets everything up for the tests of PopAnalysisTask.h/.cpp .
 */
class CoupledClusterTaskTest : public ::testing::Test {
 protected:
  CoupledClusterTaskTest() {
  }

  virtual ~CoupledClusterTaskTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
* @brief Tests restricted output.
 */
TEST_F(CoupledClusterTaskTest, minimal) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  CoupledClusterTask task(act);
  task.run();
}

} /*namespace Serenity*/
