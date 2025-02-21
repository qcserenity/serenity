/**
 * @file   CDStorageController_test.cpp
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
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
#include "integrals/CDStorageController.h"
#include "integrals/CDIntegralController.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class CDStorageControllerTest : public ::testing::Test {
 protected:
  CDStorageControllerTest() {
  }

  virtual ~CDStorageControllerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests CDStorageController.h/.cpp: Test dumping vectors on disk and load them again
 */
TEST_F(CDStorageControllerTest, testDiskDump) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  auto basCont = systemController->getBasisController();
  auto auxBasCont = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY);
  cdIntController->getACDVectors(basCont, auxBasCont);

  {
    auto cdStore = cdIntController->getStorageController("ACDAO");
    cdIntController->setDiskMode();
    cdStore->freeMemory();
    cdStore->flushFile();

    auto batchsize = cdStore->loadBatch(0);
    for (unsigned int J = 0; J < batchsize; J++) {
      auto vec = cdStore->loadVector(J);
      EXPECT_TRUE(vec);
    }
    EXPECT_EQ(batchsize, cdStore->getNVectors());
  }

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  cdIntController->cleanup();

  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
