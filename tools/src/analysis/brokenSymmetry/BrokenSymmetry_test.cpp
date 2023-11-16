/**
 * @file   BrokenSymmetry_test.cpp
 * @author Anja Massolle
 *
 * @date   17. November 2020
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
#include "analysis/brokenSymmetry/BrokenSymmetry.h"
#include "geometry/Geometry.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/BrokenSymmetryTask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/SystemAdditionTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class BrokenSymmetryTest : public ::testing::Test {
 protected:
  BrokenSymmetryTest() {
  }

  virtual ~BrokenSymmetryTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests the calculation of Magnetic Exchange Coupling constants with BS-DFT
 */
TEST_F(BrokenSymmetryTest, BSDFT) {
  auto h2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, true);
  BrokenSymmetry calcBS({h2});
  calcBS.bsDFT(1, 1);
  EXPECT_NEAR(-510.0, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-255.0, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ3UHF(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ4(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(2.0, calcBS.getS2uhfHS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2uhfBS(), 0.1);
  EXPECT_NEAR(0.13, calcBS.getSab(), 0.01);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2->getSystemPath() + "BrokenSymmetrySystem/",
                                                        "BrokenSymmetrySystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, BSDFTmethyl) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Act_def2_SVP_PBE, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_Env_def2_SVP_PBE, true);
  auto superSystemSettings = act->getSettings();
  superSystemSettings.name = "HighSpinSystem";
  superSystemSettings.charge = 0;
  superSystemSettings.spin = 0;
  auto superSystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSystemSettings);

  SystemAdditionTask<UNRESTRICTED> hsSysAdd(superSystem, {act, env});
  hsSysAdd.settings.addOccupiedOrbitals = false;
  hsSysAdd.run();

  BrokenSymmetry calcBS({superSystem});
  calcBS.bsDFT(1, 1);
  EXPECT_NEAR(-54109.0, calcBS.getJ1(), 15.0);
  EXPECT_NEAR(-27055.0, calcBS.getJ2(), 15.0);
  EXPECT_NEAR(-27005.0, calcBS.getJ3UHF(), 15.0);
  EXPECT_NEAR(2.0, calcBS.getS2uhfHS(), 0.1);
  EXPECT_NEAR(0.0, calcBS.getS2uhfBS(), 0.1);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(superSystem->getSystemPath() + "BrokenSymmetrySystem/",
                                                        "BrokenSymmetrySystem");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(superSystem->getSystemPath(), "HighSpinSystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, BSDFTload) {
  auto h2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, false);
  auto h2BS = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BS_3A, false);
  EmbeddingSettings embedding;
  BrokenSymmetry calcBS({h2}, embedding, false, false, Options::ORTHOGONALIZATION_ALGORITHMS::NONE, {h2BS});
  EXPECT_NEAR(-510.0, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-255.0, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ3UHF(), 5.0);
  EXPECT_NEAR(-501.0, calcBS.getJ4(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(2.0, calcBS.getS2uhfHS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2uhfBS(), 0.1);
  EXPECT_NEAR(0.13, calcBS.getSab(), 0.01);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2->getSystemPath() + "BrokenSymmetrySystem/",
                                                        "BrokenSymmetrySystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, FAT_NAKE) {
  auto h2act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A, true);
  auto h2env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A, true);
  EmbeddingSettings embedding;
  embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  FreezeAndThawTaskSettings fatSettings;
  BrokenSymmetry calcBS({h2act, h2env}, embedding, false, false);
  calcBS.bsFAT(fatSettings);
  EXPECT_NEAR(-348, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-174, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-342, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2act->getSystemPath() + "brokenSymmetryActSystem/",
                                                        "brokenSymmetryActSystem");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2env->getSystemPath() + "brokenSymmetryEnvSystem/",
                                                        "brokenSymmetryEnvSystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, FDE_TsOrtho) {
  auto h2act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A, false);
  auto h2env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A, false);
  auto h2BSact = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSAct_3A, false);
  auto h2BSenv = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSEnv_3A, false);
  EmbeddingSettings embedding;
  embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  BrokenSymmetry calcBS({h2act, h2env}, embedding, true, false, Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN,
                        {h2BSact, h2BSenv});
  EXPECT_NEAR(-1241, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-621, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-1223, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2act->getSystemPath() + "TestSystem_H2_6_311G_Act_3AOrtho/",
                                                        "TestSystem_H2_6_311G_Act_3AOrtho");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(
      h2BSact->getSystemPath() + "TestSystem_H2_6_311G_BSAct_3AOrtho/", "TestSystem_H2_6_311G_BSAct_3AOrtho");
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, FDE_AllOrtho) {
  auto h2act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A, false);
  auto h2env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A, false);
  auto h2BSact = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSAct_3A, false);
  auto h2BSenv = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_BSEnv_3A, false);
  EmbeddingSettings embedding;
  embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  BrokenSymmetry calcBS({h2act, h2env}, embedding, false, true, Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN,
                        {h2BSact, h2BSenv});
  EXPECT_NEAR(-417, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-208, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-410, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2BSact->getSystemPath() + "superBS/", "superBS");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2act->getSystemPath() + "superHS/", "superHS");
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, FDE_NAKE) {
  auto h2act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A, true);
  auto h2env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A, true);
  EmbeddingSettings embedding;
  embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  BrokenSymmetry calcBS({h2act, h2env}, embedding);
  calcBS.bsFDE();
  EXPECT_NEAR(-347, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-174, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-342, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2act->getSystemPath() + "brokenSymmetryActSystem/",
                                                        "brokenSymmetryActSystem");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2env->getSystemPath() + "brokenSymmetryEnvSystem/",
                                                        "brokenSymmetryEnvSystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(BrokenSymmetryTest, ISOLATED_NAKE) {
  auto h2act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Act_3A, true);
  auto h2env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_Env_3A, true);
  EmbeddingSettings embedding;
  embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::BLYP;
  BrokenSymmetry calcBS({h2act, h2env}, embedding);
  calcBS.bsIsolated();
  EXPECT_NEAR(-347, calcBS.getJ1(), 5.0);
  EXPECT_NEAR(-173, calcBS.getJ2(), 5.0);
  EXPECT_NEAR(-342, calcBS.getJ3(), 5.0);
  EXPECT_NEAR(2.0, calcBS.getS2HS(), 0.1);
  EXPECT_NEAR(1.0, calcBS.getS2BS(), 0.1);
  // Remove dummy systems created by BS task
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2act->getSystemPath() + "brokenSymmetryActSystem/",
                                                        "brokenSymmetryActSystem");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(h2env->getSystemPath() + "brokenSymmetryEnvSystem/",
                                                        "brokenSymmetryEnvSystem");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/