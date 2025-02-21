/**
 * @file LocalMP2InteractionCalculator_test.cpp
 *
 * @date 29 Apr 2020
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
#include "postHF/MPn/LocalMP2InteractionCalculator.h"           //To be tested.
#include "geometry/Geometry.h"                                  //Supersystem construction.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //LocalCorrelation settings.
#include "settings/Settings.h"                                  //Change settings of the test system.
#include "system/SystemController.h"                            //getSettings().
#include "tasks/LocalizationTask.h"                             //Orbital localization before splitting.
#include "tasks/ScfTask.h"                                      //Clean electronic structures.
#include "tasks/SystemAdditionTask.h"                           //System addition.
#include "tasks/SystemSplittingTask.h"                          //System splitting.
#include "testsupply/SystemController__TEST_SUPPLY.h"           //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class LocalMP2InteractionCalculatorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(LocalMP2InteractionCalculatorTest, WaterDimer) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  Settings supersettings = sys1->getSettings();
  supersettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  sys1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, supersettings);
  supersettings.name = "TestSystem_WaterMonTwo_Def2_SVP";
  auto sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  Settings sys2Settings = sys2->getSettings();
  sys2Settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  sys2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, sys2Settings);
  supersettings.name = "TMP_SUPERSYSTEM";
  auto supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), supersettings);
  SystemAdditionTask<RESTRICTED> add(supersystem, {sys1, sys2});
  add.settings.addOccupiedOrbitals = false;
  add.run();
  ScfTask<scfMode> scfTask(supersystem);
  scfTask.run();
  LocalizationTask loc(supersystem);
  loc.settings.splitValenceAndCore = true;
  loc.run();
  SystemSplittingTask<scfMode> split(supersystem, {sys1, sys2});
  split.run();
  Eigen::VectorXi assignment = split.getFinalAssignment();
  std::vector<unsigned int> environmentIndices;
  for (unsigned int iOcc = 0; iOcc < assignment.size(); ++iOcc)
    if (assignment(iOcc) == 1)
      environmentIndices.push_back(iOcc);
  LocalCorrelationSettings lcSettings;
  lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  LocalMP2InteractionCalculator mp2InteractionCalculator(sys1, {sys2}, lcSettings, 100, 1e-5, false, supersystem,
                                                         environmentIndices);
  double interactionEnergy = mp2InteractionCalculator.getInteractionEnergy();
  double environmentEnergy = mp2InteractionCalculator.getEnvironmentEnergy();

  EXPECT_NEAR(-0.0018756709779051543, interactionEnergy, 1e-6);
  EXPECT_NEAR(-0.20332318167847344, environmentEnergy, 1e-6);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(supersystem);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
