/**
 * @file OrthogonalizationTask_test.cpp
 *
 * @date Nov 19, 2020
 * @author Anja Massolle
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
#include "tasks/OrthogonalizationTask.h"
#include "geometry/Geometry.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class OrthogonalizationTaskTest : public ::testing::Test {
 protected:
  OrthogonalizationTaskTest() {
  }

  virtual ~OrthogonalizationTaskTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(OrthogonalizationTaskTest, loewdin_RESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT);

  OrthogonalizationTask<RESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(OrthogonalizationTaskTest, pipek_RESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env);

  OrthogonalizationTask<RESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::PIPEK;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(OrthogonalizationTaskTest, broer_RESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Act);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86_Env);

  OrthogonalizationTask<RESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::BROER;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}
TEST_F(OrthogonalizationTaskTest, loewdin_UNRESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE);

  OrthogonalizationTask<UNRESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(OrthogonalizationTaskTest, pipek_UNRESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE);

  OrthogonalizationTask<UNRESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::PIPEK;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(OrthogonalizationTaskTest, broer_UNRESTRICTED) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE);

  OrthogonalizationTask<UNRESTRICTED> orthoTask({act, env});
  orthoTask.settings.orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::BROER;
  orthoTask.run();
  ASSERT_EQ(orthoTask.isOrtho(), true);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act->getSystemPath() + (act->getSettings()).name + "Ortho/",
                                                        (act->getSettings()).name + "Ortho");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/