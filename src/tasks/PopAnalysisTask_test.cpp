/**
 * @file PopAnalysisTask_test.cpp
 *
 * @date Aug 22, 2017
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
#include "tasks/PopAnalysisTask.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class PopAnalysisTaskTest
 * @brief Sets everything up for the tests of PopAnalysisTask.h/.cpp .
 */
class PopAnalysisTaskTest : public ::testing::Test {
 protected:
  PopAnalysisTaskTest() {
  }

  virtual ~PopAnalysisTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests restricted output.
 */
TEST_F(PopAnalysisTaskTest, Mulliken_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  PopAnalysisTask<RESTRICTED> task(act);
  task.settings.algorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN;
  task.run();
}
/**
 * @test
 * @brief Tests restricted output for Hirshfeld analysis.
 */
TEST_F(PopAnalysisTaskTest, Hirshfeld_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  PopAnalysisTask<RESTRICTED> task(act);
  task.settings.algorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD;
  task.run();
  std::string hFree = act->getSettings().path + "/H_FREE/H_FREE";
  std::remove((hFree + ".settings").c_str());
  std::remove((hFree + ".xyz").c_str());
  std::remove((hFree + ".energies.unres").c_str());
  std::remove((hFree + ".orbs.unres.h5").c_str());
  std::remove((hFree + ".dmat.unres.h5").c_str());
  std::remove((act->getSettings().path + "/H_FREE").c_str());
}
/**
 * @test
 * @brief Tests restricted output for IAO analysis.
 */
TEST_F(PopAnalysisTaskTest, IAO_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  PopAnalysisTask<RESTRICTED> task(act);
  task.settings.algorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAO;
  task.run();
}
/**
 * @test
 * @brief Tests restricted output for IAO-Shell analysis.
 */
TEST_F(PopAnalysisTaskTest, IAOShell_R_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, true);
  PopAnalysisTask<RESTRICTED> task(act);
  task.settings.algorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell;
  task.run();
  std::string hFree = act->getSettings().path + "/H_FREE/H_FREE";
  std::remove((hFree + ".settings").c_str());
  std::remove((hFree + ".xyz").c_str());
  std::remove((hFree + ".energies.unres").c_str());
  std::remove((hFree + ".orbs.unres.h5").c_str());
  std::remove((hFree + ".dmat.unres.h5").c_str());
  std::remove((act->getSettings().path + "/H_FREE").c_str());
}
/**
 * @test
 * @brief Tests unrestricted output.
 */
TEST_F(PopAnalysisTaskTest, Mulliken_U_Output) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  act->setSCFMode(UNRESTRICTED);
  PopAnalysisTask<UNRESTRICTED> task(act);
  task.settings.algorithm = Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN;
  task.run();
}

} /*namespace Serenity*/
