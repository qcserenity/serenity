/**
 * @file Kernel_test.cpp
 *
 * @date Nov 09, 2017
 * @author Michael Boeckers
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
#include "postHF/LRSCF/Kernel/Kernel.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/GridController__TEST_SUPPLY.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class KernelTest
 * @brief Sets everything up for the tests of Kernel.h/.cpp .
 */
class KernelTest : public ::testing::Test {
 protected:
  KernelTest() {
  }

  virtual ~KernelTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(KernelTest, KERNEL) {
  // Actual values are checked in KernelSigmaVector_test.cpp.

  // Setup test systems
  std::vector<std::shared_ptr<SystemController>> systems(2);
  systems[0] = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE);
  systems[1] = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT);
  LRSCFTaskSettings settings;

  settings.func = CompositeFunctionals::XCFUNCTIONALS::PW91;
  settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;

  // Setup kernel
  auto kernel = std::make_shared<Kernel<Options::SCF_MODES::RESTRICTED>>(systems, systems, settings, systems.size());
  EXPECT_NO_FATAL_FAILURE(auto pp = kernel->getPP(0, 1, systems[0]->getSettings().grid.blocksize, 0);
                          auto pg = kernel->getPG(0, 1, systems[0]->getSettings().grid.blocksize, 0);
                          auto gg = kernel->getGG(0, 1, systems[0]->getSettings().grid.blocksize, 0););
}

} // namespace Serenity
