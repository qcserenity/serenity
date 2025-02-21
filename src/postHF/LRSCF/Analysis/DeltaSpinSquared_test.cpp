/**
 * @file DeltaSpinSquared_test.cpp
 *
 * @date 29 Jan. 2021
 * @author Johannes Toelle
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
#include "postHF/LRSCF/Analysis/DeltaSpinSquared.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class DeltaSpinSquaredTest
 * @brief Sets everything up for the tests of DeltaSpinSquared.h/.cpp .
 */
class DeltaSpinSquaredTest : public ::testing::Test {
 protected:
  DeltaSpinSquaredTest() {
  }

  virtual ~DeltaSpinSquaredTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(DeltaSpinSquaredTest, TDA) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  Settings settings = sys->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, settings);

  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  lrscf.settings.nEigen = 3;
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.run();
  DeltaSpinSquared<Options::SCF_MODES::UNRESTRICTED> ssquared(lrscf.getLRSCFControllers(), 3, Options::LR_METHOD::TDA,
                                                              Options::LRSCF_TYPE::ISOLATED);
  ssquared.print();
  Eigen::VectorXd ref(3);
  // values from Qchem 5.2
  ref << 2.0, 0.0, 2.0;
  for (unsigned int i = 0; i < ref.size(); i++) {
    EXPECT_NEAR(ssquared.calculateSpinSquared(i), ref(i), 1.0e-6);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DeltaSpinSquaredTest, TDDFT) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, true);
  Settings settings = sys->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_def2_SVP_ACTIVE_FDE, settings);

  std::vector<std::shared_ptr<SystemController>> active;
  active.push_back(sys);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf(active);
  lrscf.settings.nEigen = 4;
  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.run();
  DeltaSpinSquared<Options::SCF_MODES::UNRESTRICTED> ssquared(lrscf.getLRSCFControllers(), 4, Options::LR_METHOD::TDDFT,
                                                              Options::LRSCF_TYPE::ISOLATED);
  ssquared.print();
  Eigen::VectorXd ref(4);
  ref << 2.0, 0.0, 2.0, 0.0;
  for (unsigned int i = 0; i < ref.size(); i++) {
    EXPECT_NEAR(ssquared.calculateSpinSquared(i), ref(i), 1.0e-6);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */