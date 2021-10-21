/**
 * @file ExportGridTask_test.cpp
 *
 * @date Mar 6, 2017
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
#include "tasks/ExportGridTask.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <fstream>
#include <string>

namespace Serenity {

/**
 * @class ExportGridTaskTest
 * @brief Sets everything up for the tests of ExportGridTask.h/.cpp .
 */
class ExportGridTaskTest : public ::testing::Test {
 protected:
  ExportGridTaskTest() {
  }

  virtual ~ExportGridTaskTest() = default;

  inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    if (f.good()) {
      f.close();
      return true;
    }
    else {
      f.close();
      return false;
    }
  }
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests ExportGridTask.h/.cpp: Generate file.
 */
TEST_F(ExportGridTaskTest, FileGenerated_noAtomInfo) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto task = ExportGridTask(systemController);
  task.settings.withAtomInfo = false;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS.grid").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS.grid").c_str()));
}

/**
 * @test
 * @brief Tests ExportGridTask.h/.cpp: Generate file.
 */
TEST_F(ExportGridTaskTest, FileGenerated_withAtomInfo) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto task = ExportGridTask(systemController);
  task.settings.withAtomInfo = true;
  task.run();
  EXPECT_TRUE(fileExists((systemController->getSystemPath() + "TestSystem_H2_MINBAS.grid").c_str()));
  EXPECT_EQ(0, std::remove((systemController->getSystemPath() + "TestSystem_H2_MINBAS.grid").c_str()));
}

} /* namespace Serenity */
