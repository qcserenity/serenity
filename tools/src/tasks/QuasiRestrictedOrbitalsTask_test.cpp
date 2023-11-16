/**
 * @file QuasiRestrictedOrbitalsTask_test.cpp
 *
 * @date Mai 10, 2021
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
#include "tasks/QuasiRestrictedOrbitalsTask.h"        //To be tested.
#include "data/ElectronicStructure.h"                 //getEnergy();
#include "settings/Settings.h"                        //Change to HF.
#include "system/SystemController.h"                  //GetElectronicStructure.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class QuasiRestrictedOrbitalsTaskTest : public ::testing::Test {
 protected:
  QuasiRestrictedOrbitalsTaskTest() {
  }
  virtual ~QuasiRestrictedOrbitalsTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(QuasiRestrictedOrbitalsTaskTest, methyl) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  QuasiRestrictedOrbitalsTask<UNRESTRICTED> qroTask(system, {});
  qroTask.run();
  Eigen::VectorXd occupations = qroTask.getOccupationNumbers();
  EXPECT_NEAR(2.0, occupations(0), 1e-5);
  EXPECT_NEAR(2.0, occupations(1), 1e-3);
  EXPECT_NEAR(2.0, occupations(2), 1e-3);
  EXPECT_NEAR(2.0, occupations(3), 1e-2);
  EXPECT_NEAR(1.0, occupations(4), 1e-5);
  EXPECT_NEAR(0.0, occupations(5), 1e-2);
  EXPECT_NEAR(0.0, occupations(6), 1e-3);
  EXPECT_NEAR(0.0, occupations(7), 1e-3);
  EXPECT_NEAR(0.0, occupations(8), 1e-6);
  EXPECT_NEAR(0.0, occupations(9), 1e-9);
  EXPECT_NEAR(system->getElectronicStructure<UNRESTRICTED>()->getEnergy(), -39.528849014, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(QuasiRestrictedOrbitalsTaskTest, methyl_nonCanonical) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, true);
  Settings settings = system->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE, settings, 0, 1);
  QuasiRestrictedOrbitalsTask<UNRESTRICTED> qroTask(system, {});
  qroTask.settings.canonicalize = false;
  qroTask.run();
  Eigen::VectorXd occupations = qroTask.getOccupationNumbers();
  EXPECT_NEAR(2.0, occupations(0), 1e-5);
  EXPECT_NEAR(2.0, occupations(1), 1e-3);
  EXPECT_NEAR(2.0, occupations(2), 1e-3);
  EXPECT_NEAR(2.0, occupations(3), 1e-2);
  EXPECT_NEAR(1.0, occupations(4), 1e-5);
  EXPECT_NEAR(0.0, occupations(5), 1e-2);
  EXPECT_NEAR(0.0, occupations(6), 1e-3);
  EXPECT_NEAR(0.0, occupations(7), 1e-3);
  EXPECT_NEAR(0.0, occupations(8), 1e-6);
  EXPECT_NEAR(0.0, occupations(9), 1e-9);
  EXPECT_NEAR(system->getElectronicStructure<UNRESTRICTED>()->getEnergy(), -39.528849014, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(QuasiRestrictedOrbitalsTaskTest, water) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  QuasiRestrictedOrbitalsTask<UNRESTRICTED> qroTask(system, {});
  qroTask.run();
  Eigen::VectorXd occupations = qroTask.getOccupationNumbers();
  EXPECT_NEAR(2.0, occupations(0), 1e-6);
  EXPECT_NEAR(2.0, occupations(1), 1e-6);
  EXPECT_NEAR(2.0, occupations(2), 1e-6);
  EXPECT_NEAR(2.0, occupations(3), 1e-6);
  EXPECT_NEAR(2.0, occupations(4), 1e-6);
  EXPECT_NEAR(0.0, occupations(5), 1e-6);
  EXPECT_NEAR(0.0, occupations(6), 1e-6);
  EXPECT_NEAR(0.0, occupations(7), 1e-6);
  EXPECT_NEAR(0.0, occupations(8), 1e-6);
  EXPECT_NEAR(0.0, occupations(9), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(QuasiRestrictedOrbitalsTaskTest, restricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  QuasiRestrictedOrbitalsTask<RESTRICTED> qroTask(system, {});
  qroTask.run();
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/