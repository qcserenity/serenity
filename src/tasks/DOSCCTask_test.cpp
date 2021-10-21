/**
 * @file DOSCCTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 9, 2021
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
#include "tasks/DOSCCTask.h"                          //To be tested.
#include "math/FloatMaths.h"                          //Precision.
#include "settings/Settings.h"                        //Change system settings.
#include "system/SystemController.h"                  //Access to system settings.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class DOSCCTaskTest : public ::testing::Test {
 protected:
  DOSCCTaskTest() = default;

  virtual ~DOSCCTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
  void cleanUp(std::shared_ptr<SystemController> system) {
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  }
};

TEST_F(DOSCCTaskTest, ethane_bond_streching) {
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings aSettings = a->getSettings();
  aSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  aSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, aSettings);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  Settings bSettings = b->getSettings();
  bSettings.name = "b";
  bSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, bSettings);

  DOSCCTask dosccTask({a, b});
  dosccTask.settings.dosSettings = Options::DOS_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.run();
  const auto relativeEnergies = dosccTask.getRelativeEnergies();
  const auto totalEnergies = dosccTask.getTotalEnergies();

  EXPECT_NEAR(relativeEnergies(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergies(1), 0.062922261156970194, NORMAL_D);
  EXPECT_NEAR(totalEnergies(0), -79.52253707479764, NORMAL_D);
  EXPECT_NEAR(totalEnergies(1), -79.45961481364067, NORMAL_D);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(a);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(b);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
