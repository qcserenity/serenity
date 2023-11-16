/**
 * @file BasisSetTruncationTask_test.cpp
 *
 * @date Nov 15, 2018
 * @author Moritz Bensberg
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
#include "tasks/BasisSetTruncationTask.h"
#include "basis/BasisController.h"
#include "system/SystemController.h" //System definition.
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class BasisSetTruncationTaskTest BasisSetTruncationTask_test.cpp
 * @brief Sets up the test environment for the BasisSetTruncationTask.
 */
class BasisSetTruncationTaskTest : public ::testing::Test {
 protected:
  BasisSetTruncationTaskTest() = default;

  virtual ~BasisSetTruncationTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
}; // class

/**
 * @test
 * @brief Tests the net-population truncation scheme for an ethane molecule (Act:CC, Env: all H).
 *        Unrestricted tests are not necessary. The unrestricted parts of the task are already tested by different tests.
 */
TEST_F(BasisSetTruncationTaskTest, EthaneCCBond_netPop) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  BasisSetTruncationTask<SPIN> task(act);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  task.settings.netThreshold = 1e-3;
  task.run();
  EXPECT_EQ((unsigned int)34, act->getBasisController()->getNBasisFunctions());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the scaled version of the net-population truncation scheme for an ethane molecule (Act:CC, Env: all H).
 *        Unrestricted tests are not necessary. The unrestricted parts of the task are already tested by different tests.
 */
TEST_F(BasisSetTruncationTaskTest, EthaneCCBond_netPopPrim) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  BasisSetTruncationTask<SPIN> task(act);
  task.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION;
  task.settings.truncationFactor = 0.7;
  task.run();
  EXPECT_EQ((unsigned int)37, act->getBasisController()->getNBasisFunctions());
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
