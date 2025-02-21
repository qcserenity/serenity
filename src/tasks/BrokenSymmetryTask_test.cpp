/**
 * @file BrokenSymmetryTask_test.cpp
 *
 * @date Feb. 28, 2023
 * @author: Moritz Bensberg
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
#include "tasks/BrokenSymmetryTask.h" //To be tested.
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class BrokenSymmetryTaskTest : public ::testing::Test {
 protected:
  BrokenSymmetryTaskTest() = default;

  virtual ~BrokenSymmetryTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(BrokenSymmetryTaskTest, BSDFT) {
  auto h2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, true);
  BrokenSymmetryTask bsTask({h2});
  bsTask.run();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2->getSystemPath() + "BrokenSymmetrySystem/",
                                                        "BrokenSymmetrySystem");
}

// ToDo: These tests give only NaN.
// TEST_F(BrokenSymmetryTaskTest, BSDFTFDE) {
//  auto h2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, false);
//  auto h2BS = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BS_3A, false);
//  BrokenSymmetryTask bsTask({h2, h2BS});
//  bsTask.settings.embeddingScheme = Options::EMBEDDING_SCHEME::FDE;
//  bsTask.run();
//  SystemController__TEST_SUPPLY::cleanUp();
//}
//
// TEST_F(BrokenSymmetryTaskTest, BSDFTFAT) {
//  auto h2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, false);
//  auto h2BS = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BS_3A, false);
//  BrokenSymmetryTask bsTask({h2, h2BS});
//  bsTask.settings.embeddingScheme = Options::EMBEDDING_SCHEME::FAT;
//  bsTask.run();
//  SystemController__TEST_SUPPLY::cleanUp();
//}

} /* namespace Serenity */
