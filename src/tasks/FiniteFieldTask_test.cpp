/**
 * @file FiniteFieldTask_test.cpp
 *
 * @date Apr. 23, 2020
 * @author Patrick Eschenbach
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
#include "tasks/FiniteFieldTask.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class FiniteFieldTaskTest
 * @brief Sets everything up for the tests of FiniteFieldTask.h/.cpp .
 */
class FiniteFieldTaskTest : public ::testing::Test {
 protected:
  FiniteFieldTaskTest() {
  }

  virtual ~FiniteFieldTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Numerical static polarizability (FFTask) vs analytical static polarizability (LRSCFTask)
 */
TEST_F(FiniteFieldTaskTest, numerStatPolyFFvsLRSCF) {
  auto sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
  Settings settings = sysDummy->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, settings);

  FiniteFieldTask<RESTRICTED> ff(sysDummy);
  ff.settings.frequency = 0;
  ff.run();
  LRSCFTask<RESTRICTED> lrscf({sysDummy});
  lrscf.settings.frequencies = {0};
  lrscf.settings.nEigen = 0;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  double isoPoly = (1.0 / 3.0) * (std::get<1>(lrscf.getProperties()[0])).trace();
  EXPECT_NEAR(ff.getIsotropicPolarizability(), isoPoly, 5e-4);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sysDummy);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
}

/**
 * @test
 * @brief Analytical dynamic polarizability (FFTask) vs analytical dynamic polarizability (LRSCFTask)
 */
TEST_F(FiniteFieldTaskTest, analyDynPolyFFvsLRSCF) {
  auto sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
  Settings settings = sysDummy->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, settings);
  FiniteFieldTask<RESTRICTED> ff(sysDummy);
  ff.settings.frequency = 1.0;
  ff.run();
  LRSCFTask<RESTRICTED> lrscf({sysDummy});
  lrscf.settings.frequencies = {1.0};
  lrscf.settings.nEigen = 0;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  double isoPoly = (1.0 / 3.0) * (std::get<1>(lrscf.getProperties()[0])).trace();
  EXPECT_NEAR(ff.getIsotropicPolarizability(), isoPoly, 5e-7);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sysDummy);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
}

/**
 * @test
 * @brief Analytical static polarizability (FFTask) vs analytical static polarizability (LRSCFTask)
 */
TEST_F(FiniteFieldTaskTest, analyStatPolyFFvsLRSCF) {
  auto sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
  Settings settings = sysDummy->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, settings);
  FiniteFieldTask<RESTRICTED> ff(sysDummy);
  ff.settings.frequency = 0;
  ff.run();
  LRSCFTask<RESTRICTED> lrscf({sysDummy});
  lrscf.settings.frequencies = {0};
  lrscf.settings.nEigen = 0;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  double isoPoly = (1.0 / 3.0) * (std::get<1>(lrscf.getProperties()[0])).trace();
  EXPECT_NEAR(ff.getIsotropicPolarizability(), isoPoly, 5e-4);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sysDummy);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
}

/**
 * @test
 * @brief Numerical (static, static) hyperpolarizability (FFTask) vs analytical (static, static) hyperpolarizability
 * (Turbomole), 4th june 2021.
 */
TEST_F(FiniteFieldTaskTest, staticstatic) {
  auto sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
  Settings settings = sysDummy->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  sysDummy = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, settings);
  FiniteFieldTask<RESTRICTED> ff(sysDummy);
  ff.settings.frequency = 0;
  ff.settings.hyperPolarizability = true;
  ff.settings.finiteFieldStrength = 1e-3;
  ff.run();
  auto hyper = ff.getHyperPolarizability();
  Eigen::Matrix3d hypTurboX = Eigen::Matrix3d::Zero(3, 3);
  Eigen::Matrix3d hypTurboY = Eigen::Matrix3d::Zero(3, 3);
  Eigen::Matrix3d hypTurboZ = Eigen::Matrix3d::Zero(3, 3);

  hypTurboX(0, 0) = -4.1781760185159;
  hypTurboX(0, 1) = -5.3100913826039;
  hypTurboX(0, 2) = -1.4328164619266;
  hypTurboX(1, 0) = -5.3100913830111;
  hypTurboX(1, 1) = -7.4837657850723;
  hypTurboX(1, 2) = 0.39229379488152;
  hypTurboX(2, 0) = -1.4328164617963;
  hypTurboX(2, 1) = 0.39229379489145;
  hypTurboX(2, 2) = 15.077467118910;

  hypTurboY(0, 0) = -5.3100913830111;
  hypTurboY(0, 1) = -7.4837657850723;
  hypTurboY(0, 2) = 0.39229379488152;
  hypTurboY(1, 0) = -7.4837657847634;
  hypTurboY(1, 1) = -18.678407744984;
  hypTurboY(1, 2) = 2.9335904009957;
  hypTurboY(2, 0) = 0.39229379480583;
  hypTurboY(2, 1) = 2.9335904009368;
  hypTurboY(2, 2) = 11.591946193146;

  hypTurboZ(0, 0) = -1.4328164617963;
  hypTurboZ(0, 1) = 0.39229379489145;
  hypTurboZ(0, 2) = 15.077467118910;
  hypTurboZ(1, 0) = 0.39229379480583;
  hypTurboZ(1, 1) = 2.9335904009368;
  hypTurboZ(1, 2) = 11.591946193146;
  hypTurboZ(2, 0) = 15.077467118505;
  hypTurboZ(2, 1) = 11.591946193027;
  hypTurboZ(2, 2) = -7.3691418536523;

  double accuracy = 1e-1;
  EXPECT_NEAR(hyper[0](0, 0), hypTurboX(0, 0), accuracy);
  EXPECT_NEAR(hyper[0](1, 0), hypTurboX(1, 0), accuracy);
  EXPECT_NEAR(hyper[0](2, 0), hypTurboX(2, 0), accuracy);
  EXPECT_NEAR(hyper[0](0, 1), hypTurboX(0, 1), accuracy);
  EXPECT_NEAR(hyper[0](1, 1), hypTurboX(1, 1), accuracy);
  EXPECT_NEAR(hyper[0](0, 2), hypTurboX(0, 2), accuracy);
  EXPECT_NEAR(hyper[0](2, 2), hypTurboX(2, 2), accuracy);

  EXPECT_NEAR(hyper[1](0, 0), hypTurboY(0, 0), accuracy);
  EXPECT_NEAR(hyper[1](1, 0), hypTurboY(1, 0), accuracy);
  EXPECT_NEAR(hyper[1](0, 1), hypTurboY(0, 1), accuracy);
  EXPECT_NEAR(hyper[1](1, 1), hypTurboY(1, 1), accuracy);
  EXPECT_NEAR(hyper[1](2, 1), hypTurboY(2, 1), accuracy);
  EXPECT_NEAR(hyper[1](1, 2), hypTurboY(1, 2), accuracy);
  EXPECT_NEAR(hyper[1](2, 2), hypTurboY(2, 2), accuracy);

  EXPECT_NEAR(hyper[2](0, 0), hypTurboZ(0, 0), accuracy);
  EXPECT_NEAR(hyper[2](2, 0), hypTurboZ(2, 0), accuracy);
  EXPECT_NEAR(hyper[2](1, 1), hypTurboZ(1, 1), accuracy);
  EXPECT_NEAR(hyper[2](2, 1), hypTurboZ(2, 1), accuracy);
  EXPECT_NEAR(hyper[2](0, 2), hypTurboZ(0, 2), accuracy);
  EXPECT_NEAR(hyper[2](1, 2), hypTurboZ(1, 2), accuracy);
  EXPECT_NEAR(hyper[2](2, 2), hypTurboZ(2, 2), accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sysDummy);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF);
}

} /*namespace Serenity*/
