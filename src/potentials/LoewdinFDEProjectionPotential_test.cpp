/**
 * @file LoewdinFDEProjectionPotential_test.cpp
 *
 * @date Jun 23, 2024
 * @author Denis G. Artiukhin
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
#include "potentials/LoewdinFDEProjectionPotential.h"
#include "data/ElectronicStructure.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {

/**
 * @class LoewdinFDEProjectionPotentialTest
 * @brief Sets everything up for the tests of LoewdinFDEProjectionPotential.h/.cpp .
 */
class LoewdinFDEProjectionPotentialTest : public ::testing::Test {
 protected:
  LoewdinFDEProjectionPotentialTest()
    : act(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      env(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~LoewdinFDEProjectionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// subsystems
  std::shared_ptr<SystemController> act;
  std::shared_ptr<SystemController> env;

  void remove_tmp_folders(const std::string& str) {
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res")
                    .c_str());

    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                       "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5")
                    .c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
    std::remove((str + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  }
};

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests restricted fock matrix of the 0th order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, restricted_0th_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 0;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  EXPECT_NEAR(0.0, f(0, 0), 1e-7);
  EXPECT_NEAR(0.0, f(0, 1), 1e-7);
  EXPECT_NEAR(0.0, f(1, 0), 1e-7);
  EXPECT_NEAR(0.0, f(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests restricted fock matrix of the 1st order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, restricted_1st_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 1;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  Eigen::Vector4d res(0.00642235754183966, 0.0009681368515132, 0.0009681368515132, -0.0404906949935935);
  EXPECT_NEAR(res(0), f(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests unrestricted fock matrix of the 1st order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, unrestricted_1st_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 1;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  Eigen::Vector4d res(0.00642235754183966, 0.0009681368515132, 0.0009681368515132, -0.0404906949935935);

  EXPECT_NEAR(res(0), f.alpha(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.alpha(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.alpha(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.alpha(1, 1), 1e-7);

  EXPECT_NEAR(res(0), f.beta(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.beta(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.beta(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.beta(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests restricted fock matrix of the 2nd order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, restricted_2nd_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 2;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  Eigen::Vector4d res(0.11325384669659103, 0.1832069689953866, 0.1832069689953866, 0.17186081939320921);

  EXPECT_NEAR(res(0), f(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests unrestricted fock matrix of the 2nd order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, unrestricted_2nd_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 2;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  Eigen::Vector4d res(0.11665071577131497, 0.18868078276361286, 0.18868078276361286, 0.17690252979044441);

  EXPECT_NEAR(res(0), f.alpha(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.alpha(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.alpha(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.alpha(1, 1), 1e-7);

  EXPECT_NEAR(res(0), f.beta(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.beta(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.beta(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.beta(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LoewdinFDEProjectionPotential.h/.cpp: Tests unrestricted fock matrix of the 3rd order expansion.
 */
TEST_F(LoewdinFDEProjectionPotentialTest, unrestricted_3rd_loewdin) {
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.loewdinOrder = 3;

  LoewdinFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED> loewdinFDEPotential(act, {env}, settings);
  auto f = loewdinFDEPotential.getMatrix();
  Eigen::Vector4d res(0.11507109860572058, 0.18019665462149637, 0.18019665462149637, 0.1428909120814818);

  EXPECT_NEAR(res(0), f.alpha(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.alpha(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.alpha(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.alpha(1, 1), 1e-7);

  EXPECT_NEAR(res(0), f.beta(0, 0), 1e-7);
  EXPECT_NEAR(res(1), f.beta(0, 1), 1e-7);
  EXPECT_NEAR(res(2), f.beta(1, 0), 1e-7);
  EXPECT_NEAR(res(3), f.beta(1, 1), 1e-7);

  remove_tmp_folders(env->getSystemPath());
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
