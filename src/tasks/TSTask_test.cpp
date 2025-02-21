/**
 * @file TSTask_test.cpp
 * @author: Jan Unsleber
 *
 * @date Sep 20, 2017
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
#include "tasks/TSTask.h"
#include "geometry/Geometry.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class TSTaskTest : public ::testing::Test {
 protected:
  TSTaskTest()
    : _system(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS)) {
  }
  virtual ~TSTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /* The systems. */
  std::shared_ptr<SystemController> _system;
};

/**
 * @test TSTaskTest
 * @brief Optimizes water into linear geometry using stored hessian guess.
 */
TEST_F(TSTaskTest, LinearWater) {
  /*
   * Copy hessian file
   */
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/";
  }
  else {
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set." << std::endl;
    assert(false);
  }
  std::ifstream src(
      (pathToTestsResources + "TestSystem_WATER_DISTORTED_MINBAS/TestSystem_WATER_DISTORTED_MINBAS.hess.h5").c_str(),
      std::ios::binary);
  std::ofstream dest((_system->getSystemPath() + "TestSystem_WATER_DISTORTED_MINBAS.hess.h5").c_str(), std::ios::binary);
  dest << src.rdbuf();
  dest.flush();

  /*
   * actual test
   */
  TSTask tssearch(_system, {});
  tssearch.run();
  EXPECT_TRUE(_system->getGeometry()->isLinear());
  EXPECT_EQ(0, std::remove((_system->getSystemPath() + "TestSystem_WATER_DISTORTED_MINBAS.hess.h5").c_str()));
  EXPECT_EQ(0, std::remove((_system->getSystemPath() + "TestSystem_WATER_DISTORTED_MINBAS.trj").c_str()));
}

} /* namespace Serenity */
