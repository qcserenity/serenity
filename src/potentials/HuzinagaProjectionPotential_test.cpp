/**
 * @file HuzinagaProjectionPotential_test.cpp
 *
 * @date Feb 23, 2018
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
#include "potentials/HuzinagaProjectionPotential.h"
#include "geometry/Geometry.h"
#include "potentials/NAddFuncPotential.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {

/**
 * @class HuzinagaProjectionPotentialTest
 * @brief Sets everything up for the tests of HuzinagaProjectionPotential.h/.cpp .
 */
class HuzinagaProjectionPotentialTest : public ::testing::Test {
 protected:
  HuzinagaProjectionPotentialTest() {
  }

  virtual ~HuzinagaProjectionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests HuzinagaProjectionPotential.h/.cpp: Tests restricted fock matrix.
 */
TEST_F(HuzinagaProjectionPotentialTest, restricted_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  HuzinagaProjectionPotential<Options::SCF_MODES::RESTRICTED> huzFDEPotential(act, {env}, settings);
  auto f = huzFDEPotential.getMatrix();
  EXPECT_NEAR(0.052226824418551698, f(0, 0), 1e-7);
  EXPECT_NEAR(0.11198231040394366, f(0, 1), 1e-7);
  EXPECT_NEAR(0.11198231040394366, f(1, 0), 1e-7);
  EXPECT_NEAR(0.22587526848483189, f(1, 1), 1e-7);

  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaProjectionPotential.h/.cpp: Tests unrestricted fock matrix.
 */
TEST_F(HuzinagaProjectionPotentialTest, unrestricted_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  HuzinagaProjectionPotential<Options::SCF_MODES::UNRESTRICTED> huzFDEPotential(act, {env}, settings);
  auto f = huzFDEPotential.getMatrix();
  EXPECT_NEAR(0.05166910663948323, f.alpha(0, 0), 1e-7);
  EXPECT_NEAR(0.10986760211800925, f.alpha(0, 1), 1e-7);
  EXPECT_NEAR(0.10986760211800925, f.alpha(1, 0), 1e-7);
  EXPECT_NEAR(0.21856343934609146, f.alpha(1, 1), 1e-7);

  EXPECT_NEAR(0.05166910663948323, f.beta(0, 0), 1e-7);
  EXPECT_NEAR(0.10986760211800925, f.beta(0, 1), 1e-7);
  EXPECT_NEAR(0.10986760211800925, f.beta(1, 0), 1e-7);
  EXPECT_NEAR(0.21856343934609146, f.beta(1, 1), 1e-7);
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaProjectionPotential.h/.cpp: Tests restricted fock matrix; with completely truncated projector.
 */
TEST_F(HuzinagaProjectionPotentialTest, restricted_truncProjection) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.truncateProjector = true;
  settings.projecTruncThresh = 100;
  HuzinagaProjectionPotential<SPIN> huzFDEPotential(act, {env}, settings);

  auto f = huzFDEPotential.getMatrix();
  EXPECT_NEAR(0.0, f(0, 0), 1e-7);
  EXPECT_NEAR(0.0, f(0, 1), 1e-7);
  EXPECT_NEAR(0.0, f(1, 0), 1e-7);
  EXPECT_NEAR(0.0, f(1, 1), 1e-7);

  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
