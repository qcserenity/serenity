/**
 * @file Besley_test.cpp
 *
 * @date Dec 21, 2018
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Tools/Besley.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class BesleyTest
 * @brief Sets everything up for the tests of Besley.h/.cpp .
 */
  class BesleyTest : public ::testing::Test {
  protected:
    BesleyTest() {
    }

    virtual ~BesleyTest() = default;

    static void TearDownTestCase() {
      SystemController__TEST_SUPPLY::cleanUp();
    }
  };

/**
 * @test
 * @brief Tests Besley.h/.cpp.
 */
  TEST_F(BesleyTest, res) {
    const auto SCFMode = Options::SCF_MODES::RESTRICTED;
    auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, true);
    //Trigger SCF
    systemController->getElectronicStructure<SCFMode>();

    //Define settings
    LRSCFTaskSettings settings;
    settings.besleyAtoms = 1;
    settings.besleyCutoff = {0.25, 0.6};

    Besley<SCFMode> besley(systemController, settings.besleyAtoms, settings.besleyCutoff);
    auto besleyList = besley.getWhiteList();

    for_spin(besleyList){
      EXPECT_EQ(besleyList_spin.size(), (unsigned int) 6);
      EXPECT_EQ(besleyList[0], (unsigned int) 1);
      EXPECT_EQ(besleyList[1], (unsigned int) 4);
      EXPECT_EQ(besleyList[2], (unsigned int) 5);
      EXPECT_EQ(besleyList[3], (unsigned int) 6);
      EXPECT_EQ(besleyList[4], (unsigned int) 7);
      EXPECT_EQ(besleyList[5], (unsigned int) 8);
    };
  }

} /* namespace Serenity */