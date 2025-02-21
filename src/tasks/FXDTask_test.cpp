/**
 * @file FXDTask_test.cpp
 *
 * @date Apr 29, 2019
 * @author J. Toelle
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
#include "tasks/FXDTask.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <iomanip>

namespace Serenity {

/**
 * @class LRSCFTask_test
 * @brief Sets everything up for the tests of LRSCFTask.h/.cpp .
 */
class FXDTaskTest : public ::testing::Test {
 protected:
  FXDTaskTest() {
  }

  virtual ~FXDTaskTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(FXDTaskTest, FED) {
  auto he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, true);
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings = he2->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, settings);
  // SCF
  ScfTask<SCFMode> scf(he2);
  scf.run();
  // TDDFT
  LRSCFTask<SCFMode> lrscf({he2});
  lrscf.settings.nEigen = 3;
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  // FED
  FXDTask<SCFMode> fxd(he2);
  fxd.settings.donoratoms = {0, 0};
  fxd.settings.acceptoratoms = {1, 1};
  fxd.settings.FED = true;
  fxd.settings.loewdinpopulation = false;
  fxd.run();
  Eigen::VectorXd qchem_FED(6);
  qchem_FED << 0.0, -0.964964, 0.0, 0.0, 0.0, 0.0;
  auto fed = fxd._fed_matrix;
  // CIS results qchem 5.2.1 settings same as serenity
  unsigned int counter = 0;
  for (unsigned int n = 0; n < fed[0].cols(); n++) {
    for (unsigned int m = n; m < fed[0].cols(); m++) {
      EXPECT_NEAR(std::abs(fed[0](n, m) - fed[1](n, m)), std::abs(qchem_FED(counter)), 1.0e-6);
      counter += 1;
    }
  }
}
TEST_F(FXDTaskTest, FCD) {
  auto he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, true);
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings = he2->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, settings);
  // SCF
  ScfTask<SCFMode> scf(he2);
  scf.run();
  // TDDFT
  LRSCFTask<SCFMode> lrscf({he2});
  lrscf.settings.nEigen = 3;
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  auto lrscfContr = lrscf.getLRSCFControllers();
  // FED
  FXDTask<SCFMode> fxd(he2);
  fxd.settings.donoratoms = {0, 0};
  fxd.settings.acceptoratoms = {1, 1};
  fxd.settings.FCD = true;
  fxd.settings.loewdinpopulation = false;
  fxd.run();
  Eigen::VectorXd qchem_FCD(6);
  qchem_FCD << 0.0, -0.913656, 0.0, 0.0, 0.0, 0.0;
  auto fcd = fxd._fcd_matrix;
  // CIS results qchem 5.2.1 settings same as serenity
  unsigned int counter = 0;
  for (unsigned int n = 0; n < fcd[0].cols(); n++) {
    for (unsigned int m = n; m < fcd[0].cols(); m++) {
      EXPECT_NEAR(std::abs(fcd[0](n, m) - fcd[1](n, m)), std::abs(qchem_FCD(counter)), 1.0e-6);
      counter += 1;
    }
  }
}
TEST_F(FXDTaskTest, FXD) {
  auto he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, true);
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  Settings settings = he2->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  he2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::He2_def2SVP_BP86, settings);
  // SCF
  ScfTask<SCFMode> scf(he2);
  scf.run();
  // TDDFT
  LRSCFTask<SCFMode> lrscf({he2});
  lrscf.settings.nEigen = 3;
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.run();
  auto lrscfContr = lrscf.getLRSCFControllers();
  // FXD
  FXDTask<SCFMode> fxd(he2);
  fxd.settings.donoratoms = {0, 0};
  fxd.settings.acceptoratoms = {1, 1};
  fxd.settings.multistateFXD = true;
  fxd.settings.states = 2;
  fxd.run();
  auto fxdMat = fxd._fxd_matrix;
  // Reference from Serenity
  EXPECT_NEAR(std::abs(fxdMat(1, 0)), std::abs(0.074467), 1.0e-6);
}

TEST_F(FXDTaskTest, FXD2) {
  auto hene = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::HeNe_def2SVP_BP86, true);
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  // SCF
  ScfTask<SCFMode> scf(hene);
  scf.run();
  // TDDFT
  LRSCFTask<SCFMode> lrscf({hene});
  lrscf.settings.nEigen = 8;
  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.blockAveThreshold = 1e-16;
  lrscf.run();
  auto lrscfContr = lrscf.getLRSCFControllers();
  // FXD
  FXDTask<SCFMode> fxd(hene);
  fxd.settings.donoratoms = {0, 0};
  fxd.settings.acceptoratoms = {1, 1};
  fxd.settings.multistateFXD = true;
  fxd.settings.states = 8;
  fxd.settings.loewdinpopulation = true;
  fxd.run();
  auto fxdMat = fxd._fxd_matrix;
  // Reference from Serenity 9.6.20
  std::cout << std::fixed << std::setprecision(9);

  Eigen::VectorXd ref(8);
  ref << 1.457480519, 1.465554476, 1.479825045, 1.725460908, 1.742723550, 1.744081862, 1.746990445, 1.749071772;
  // Checks the diabatic excitation energies
  for (unsigned int i = 0; i < fxdMat.cols(); i++) {
    EXPECT_NEAR(fxdMat(i, i), ref(i), 1.0e-6);
  }
}
} // namespace Serenity
