/**
 * @file HoffmannProjectionPotential_test.cpp
 *
 * @author Moritz Bensberg
 * @date Aug 26, 2019
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
#include "potentials/HoffmannProjectionPotential.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {
class HoffmannProjectionPotentialTest : public ::testing::Test {
 protected:
  HoffmannProjectionPotentialTest() {
  }

  virtual ~HoffmannProjectionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Tests HoffmannProjectionPotential.h/.cpp: Tests restricted fock matrix with unrestricted environment.
 */
TEST_F(HoffmannProjectionPotentialTest, restricted_unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  env->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  HoffmannProjectionPotential<Options::SCF_MODES::RESTRICTED> hoffPotential(act, {env}, settings);
  auto f = hoffPotential.getMatrix();
  EXPECT_NEAR(0.041308143680642548, f(0, 0), 1e-7);
  EXPECT_NEAR(0.081439053107599857, f(0, 1), 1e-7);
  EXPECT_NEAR(0.081439053107599857, f(1, 0), 1e-7);
  EXPECT_NEAR(0.14191034391458512, f(1, 1), 1e-7);
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.settings")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.xyz")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.energies.res")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.orbs.res.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.dmat.res.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/out").c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests HoffmannProjectionPotential.h/.cpp: Tests unrestricted fock matrix.
 */
TEST_F(HoffmannProjectionPotentialTest, unrestricted) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  env->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  HoffmannProjectionPotential<Options::SCF_MODES::UNRESTRICTED> hoffPotential(act, {env}, settings);
  auto f = hoffPotential.getMatrix();
  EXPECT_NEAR(0.039352521735825291, f.alpha(0, 0), 1e-7);
  EXPECT_NEAR(0.077029478641347057, f.alpha(0, 1), 1e-7);
  EXPECT_NEAR(0.077029478641347057, f.alpha(1, 0), 1e-7);
  EXPECT_NEAR(0.13101124211176896, f.alpha(1, 1), 1e-7);

  EXPECT_NEAR(0.039352521735825291, f.beta(0, 0), 1e-7);
  EXPECT_NEAR(0.077029478641347057, f.beta(0, 1), 1e-7);
  EXPECT_NEAR(0.077029478641347057, f.beta(1, 0), 1e-7);
  EXPECT_NEAR(0.13101124211176896, f.beta(1, 1), 1e-7);
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
 * @brief Tests HoffmannProjectionPotential.h/.cpp: Tests restricted fock matrix (Hoffmann) with unrestricted
 * environment.
 */
TEST_F(HoffmannProjectionPotentialTest, restricted_longRangeKinetic) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.longRangeNaddKinFunc = CompositeFunctionals::KINFUNCTIONALS::LLP91K;
  std::vector<std::shared_ptr<Atom>> superSystemAtoms;
  superSystemAtoms.push_back(act->getGeometry()->getAtoms()[0]);
  superSystemAtoms.push_back(act->getGeometry()->getAtoms()[1]);
  superSystemAtoms.push_back(env->getGeometry()->getAtoms()[0]);
  superSystemAtoms.push_back(env->getGeometry()->getAtoms()[1]);
  auto superGeom = std::make_shared<Geometry>(superSystemAtoms);
  auto supersystemGrid =
      AtomCenteredGridControllerFactory::produce(superGeom, act->getSettings().grid, Options::GRID_PURPOSES::DEFAULT);

  HoffmannProjectionPotential<Options::SCF_MODES::RESTRICTED> hoffPotential(act, {env}, settings, nullptr, false,
                                                                            supersystemGrid);
  auto f = hoffPotential.getMatrix();
  EXPECT_NEAR(0.011251152747563315, f(0, 0), 1e-7);
  EXPECT_NEAR(0.014840439061002785, f(0, 1), 1e-7);
  EXPECT_NEAR(0.014840439061002785, f(1, 0), 1e-7);
  EXPECT_NEAR(0.051706332270953115, f(1, 1), 1e-7);
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.settings")
                  .c_str());
  std::remove((env->getSystemPath() +
               "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.xyz")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.energies.res")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.orbs.res.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/"
                                      "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.dmat.res.h5")
                  .c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/out").c_str());
  std::remove((env->getSystemPath() + "TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HoffmannProjectionPotentialTest.h/.cpp: Tests restricted fock matrix; with completely truncated
 * projector.
 */
TEST_F(HoffmannProjectionPotentialTest, restricted_truncProjection) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  EmbeddingSettings settings;
  settings.naddXCFunc = act->getSettings().dft.functional;
  settings.truncateProjector = true;
  settings.projecTruncThresh = 100;
  HoffmannProjectionPotential<SPIN> hoffPotential(act, {env}, settings);

  auto f = hoffPotential.getMatrix();
  EXPECT_NEAR(0.0, f(0, 0), 1e-7);
  EXPECT_NEAR(0.0, f(0, 1), 1e-7);
  EXPECT_NEAR(0.0, f(1, 0), 1e-7);
  EXPECT_NEAR(0.0, f(1, 1), 1e-7);

  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
