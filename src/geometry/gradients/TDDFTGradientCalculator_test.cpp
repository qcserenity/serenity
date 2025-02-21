/**
 * @file   TDDFTGradientCalculator_test.cpp
 *
 * @date   Apr 25, 2024
 * @author Anton Rikus
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
#include "geometry/gradients/TDDFTGradientCalculator.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @test
 * @brief Tests TDDFTGradientCalculator: calculation of the derivative
 *        of the energy of excited states w.r.t nuclear coordinates
 */
class TDDFTGradientsTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(TDDFTGradientsTest, CIS_RESTRICTED) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 May 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.002059035, 0.000000000, -0.070291877, -0.067343807, 0.000000000, 0.020248674, 0.065284772,
      -0.000000000, 0.050043204;
  orcaref[1] << -0.001259956, -0.000000000, -0.092410449, -0.089564511, 0.000000000, 0.022789369, 0.090824466,
      -0.000000000, 0.069621080;
  orcaref[2] << 0.058613919, -0.000000000, -0.069834231, -0.052209985, 0.000000000, 0.074742982, -0.006403934,
      0.000000000, -0.004908751;
  orcaref[3] << 0.061046672, -0.000000000, -0.095174912, -0.076048319, -0.000000000, 0.083675072, 0.015001647,
      0.000000000, 0.011499840;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-6);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, CIS_UNRESTRICTED_Water) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();
  // ORCA 5.0.3 May 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.004909112, -0.000000000, -0.082926927, -0.078804715, 0.000000000, 0.026283222, 0.073895604,
      -0.000000000, 0.056643705;
  orcaref[1] << 0.002059035, 0.000000000, -0.070291877, -0.067343807, 0.000000000, 0.020248674, 0.065284772,
      -0.000000000, 0.050043204;
  orcaref[2] << 0.051634637, -0.000000000, -0.125253148, -0.107540083, -0.000000000, 0.082400151, 0.055905446,
      0.000000000, 0.042852998;
  orcaref[3] << -0.003529172, -0.000000000, -0.106313072, -0.103579295, 0.000000000, 0.024209653, 0.107108468,
      0.000000000, 0.082103419;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-6);
      }
    }
  }
}

// an open-shell molecule
TEST_F(TDDFTGradientsTest, CIS_UNRESTRICTED_O2) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_TRIP_CIS);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDA;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.excGradList = {1, 2, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);
  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();
  // ORCA 5.0.3 May 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.502152560;
  orcaref[0](1, 2) = -0.502152560;
  orcaref[1](0, 2) = 0.502152560;
  orcaref[1](1, 2) = -0.502152560;
  orcaref[2](0, 2) = 0.508620656;
  orcaref[2](1, 2) = -0.508620656;
  orcaref[3](0, 2) = 0.351865778;
  orcaref[3](1, 2) = -0.351865778;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-6);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDA_RESTRICTED) {
  Settings settings;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.incrementalSteps = 1;
  settings.basis.integralThreshold = 1e-64;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.scf.energyThreshold = 1e-10;
  settings.scf.rmsdThreshold = 1e-9;
  settings.scf.diisThreshold = 1e-9;
  auto system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs_DFT, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr, 1e-6);
  // calculate the gradients of the third and fourth excited state
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();
  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(2, Eigen::MatrixXd::Zero(3, 3));
  orcaref[0] << -0.078881464, 0.016904923, -0.049663871, 0.078460126, -0.020793127, 0.048846000, 0.000419982,
      0.003888665, 0.000819002;
  orcaref[1] << -0.098662432, 0.041417379, -0.059299716, 0.102405596, -0.006755486, 0.066587698, -0.003744552,
      -0.034661583, -0.007286938;

  for (unsigned i = 0; i < 2; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-6);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDA_UNRESTRICTED) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDA);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDA;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.conv = 5e-8;
  lrscfSettings.excGradList = {1, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);

  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr, 1e-6);
  // calculate the gradients of the first, third and fourth excited state
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();
  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(3, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.323609976;
  orcaref[0](1, 2) = -0.323609976;
  orcaref[1](0, 2) = 0.325819447;
  orcaref[1](1, 2) = -0.325819447;
  orcaref[2](0, 2) = 0.190837169;
  orcaref[2](1, 2) = -0.190837169;

  for (unsigned i = 0; i < 2; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-5);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDHF_RESTRICTED) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.002083470, 0.000000000, -0.072665648, -0.069629754, 0.000000000, 0.020888929, 0.067546283,
      -0.000000000, 0.051776719;
  orcaref[1] << -0.001649941, -0.000000000, -0.094371793, -0.091559810, 0.000000000, 0.022922284, 0.093209751,
      -0.000000000, 0.071449509;
  orcaref[2] << 0.058631151, -0.000000000, -0.073214992, -0.055470231, -0.000000000, 0.075637878, -0.003160920,
      0.000000000, -0.002422886;
  orcaref[3] << 0.059651343, -0.000000000, -0.098518587, -0.079639657, 0.000000000, 0.083196258, 0.019988313,
      0.000000000, 0.015322329;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-6);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDHF_UNRESTRICTED) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.conv = 1e-7;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.005228602, -0.000000000, -0.086766565, -0.082429565, 0.000000000, 0.027589207, 0.077200964,
      -0.000000000, 0.059177358;
  orcaref[1] << 0.002083473, 0.000000000, -0.072665647, -0.069629751, 0.000000000, 0.020888931, 0.067546278,
      -0.000000000, 0.051776715;
  orcaref[2] << 0.049277990, -0.000000000, -0.146825561, -0.128984457, -0.000000000, 0.085728560, 0.079706467,
      0.000000000, 0.061097000;
  orcaref[3] << -0.004408408, -0.000000000, -0.110128660, -0.107492263, 0.000000000, 0.024351802, 0.111900671,
      0.000000000, 0.085776858;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-6);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDHF_UNRESTRICTED_O2) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDHF);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDDFT;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.conv = 5e-8;
  lrscfSettings.excGradList = {1, 2, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);

  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 August 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.826801598;
  orcaref[0](1, 2) = -0.826801598;
  orcaref[1](0, 2) = 0.826801598;
  orcaref[1](1, 2) = -0.826801598;
  orcaref[2](0, 2) = 0.633532848;
  orcaref[2](1, 2) = -0.633532848;
  orcaref[3](0, 2) = 0.271711699;
  orcaref[3](1, 2) = -0.271711699;

  for (unsigned i = 0; i < 3; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDDFT_RESTRICTED) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
  settings.basis.label = "6-31gs";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.000086346, 0.000000000, -0.105806416, -0.102152101, 0.000000000, 0.027569669, 0.102064870,
      -0.000000000, 0.078236758;
  orcaref[1] << -0.001763458, 0.000000000, -0.126038749, -0.122167330, -0.000000000, 0.031038605, 0.123929954,
      -0.000000000, 0.095000210;
  orcaref[2] << 0.059113357, 0.000000000, -0.078319736, -0.060274583, -0.000000000, 0.077429630, 0.001160130,
      0.000000000, 0.000889976;
  orcaref[3] << 0.058658913, 0.000000000, -0.108374753, -0.089415139, -0.000000000, 0.084797776, 0.030755213,
      -0.000000000, 0.023576936;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-6);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDDFT_UNRESTRICTED) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_TRIP_DEF2_SVP_TDDFT);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDDFT;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.grid.accuracy = 5;
  lrscfSettings.grid.smallGridAccuracy = 5;
  lrscfSettings.conv = 5e-8;
  lrscfSettings.excGradList = {1, 2, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);

  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 June 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.323814350;
  orcaref[0](1, 2) = -0.323814350;
  orcaref[1](0, 2) = 0.323814350;
  orcaref[1](1, 2) = -0.323814350;
  orcaref[2](0, 2) = 0.325819447;
  orcaref[2](1, 2) = -0.325819447;
  orcaref[3](0, 2) = 0.189981213;
  orcaref[3](1, 2) = -0.189981213;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-6);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDA_R_GGA) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "6-31g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 July 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.007097305, 0.000000000, -0.108675679, -0.103097090, -0.000000000, 0.035084169, 0.096001760,
      0.000000000, 0.073595694;
  orcaref[1] << 0.071549562, -0.000000000, -0.062000747, -0.041281816, 0.000000000, 0.085197703, -0.030266563,
      0.000000000, -0.023193661;
  orcaref[2] << 0.004891924, 0.000000000, -0.129014326, -0.123316675, -0.000000000, 0.038239497, 0.118426651,
      0.000000000, 0.090778607;
  orcaref[3] << 0.067834081, 0.000000000, -0.097230224, -0.076274144, -0.000000000, 0.090764249, 0.008441650,
      -0.000000000, 0.006469270;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDA_U_GGA) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDA);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDA;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.grid.accuracy = 5;
  lrscfSettings.grid.smallGridAccuracy = 5;
  lrscfSettings.conv = 5e-8;
  lrscfSettings.excGradList = {1, 2, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);

  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 August 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.421639122;
  orcaref[0](1, 2) = -0.421639122;
  orcaref[1](0, 2) = 0.421639091;
  orcaref[1](1, 2) = -0.421639091;
  orcaref[2](0, 2) = 0.425229981;
  orcaref[2](1, 2) = -0.425229981;
  orcaref[3](0, 2) = 0.281782382;
  orcaref[3](1, 2) = -0.281782382;

  for (unsigned i = 0; i < 3; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDDFT_R_GGA) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "6-31g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 July 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.007412906, 0.000000000, -0.109018086, -0.103346989, -0.000000000, 0.035478105, 0.095936092,
      0.000000000, 0.073544210;
  orcaref[1] << 0.071016553, 0.000000000, -0.066935043, -0.046186967, 0.000000000, 0.085965218, -0.024828302,
      -0.000000000, -0.019026804;
  orcaref[2] << 0.004875872, -0.000000000, -0.129165882, -0.123468880, -0.000000000, 0.038264291, 0.118594915,
      -0.000000000, 0.090905381;
  orcaref[3] << 0.067336420, 0.000000000, -0.098522848, -0.077652245, -0.000000000, 0.090619571, 0.010317391,
      -0.000000000, 0.007906571;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDDFT_U_GGA) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_TRIP_6_31G_PBE_TDDFT);

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.method = Options::LR_METHOD::TDDFT;
  lrscfSettings.densFitJ = Options::DENS_FITS::NONE;
  lrscfSettings.nEigen = 4;
  lrscfSettings.grid.accuracy = 5;
  lrscfSettings.grid.smallGridAccuracy = 5;
  lrscfSettings.conv = 5e-8;
  lrscfSettings.excGradList = {1, 2, 3, 4};
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscfSettings);

  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 August 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.422137340;
  orcaref[0](1, 2) = -0.422137340;
  orcaref[1](0, 2) = 0.422136834;
  orcaref[1](1, 2) = -0.422136834;
  orcaref[2](0, 2) = 0.425229981;
  orcaref[2](1, 2) = -0.425229981;
  orcaref[3](0, 2) = 0.282329416;
  orcaref[3](1, 2) = -0.282329416;

  for (unsigned i = 0; i < 3; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
}

TEST_F(TDDFTGradientsTest, TDA_UvsR_GGA) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "6-31g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.scf.energyThreshold = 1e-10;
  settings.scf.rmsdThreshold = 1e-9;
  settings.scf.diisThreshold = 1e-9;
  auto rsystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({rsystem});
  auto rlrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(rsystem, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 3;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(rlrscfContr, 1e-6);

  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  auto usystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> ulrscf({usystem});
  auto ulrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(usystem, ulrscf.settings);

  ulrscf.settings.method = Options::LR_METHOD::TDA;
  ulrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  ulrscf.settings.nEigen = 3;
  ulrscf.settings.grid.accuracy = 7;
  ulrscf.settings.grid.smallGridAccuracy = 7;
  ulrscf.settings.conv = 5e-8;
  ulrscf.run();
  ulrscf.settings.excGradList = {2};
  TDDFTGradientCalculator<UNRESTRICTED> utddftgrads(ulrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> rexcgrads = tddftgrads.calculateGradients();
  std::vector<Eigen::MatrixXd> uexcgrads = utddftgrads.calculateGradients();

  for (unsigned i = 0; i < 2; i++) {
    for (unsigned j = 0; j < 3; j++) {
      EXPECT_NEAR(rexcgrads[0](i, j), uexcgrads[0](i, j), 5e-6);
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(rsystem);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(usystem);
}

TEST_F(TDDFTGradientsTest, TDDFT_R_GGA_RangeSeparatedHybrid) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.label = "6-31g";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto system =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::NONE;
  lrscf.settings.nEigen = 4;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3, 4};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr, 1e-6);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 5.0.3 October 2024
  std::vector<Eigen::MatrixXd> orcaref(4, Eigen::MatrixXd(3, 3));
  orcaref[0] << 0.009021213, -0.000000000, -0.100670595, -0.094873187, -0.000000000, 0.034864025, 0.085857044,
      -0.000000000, 0.065811556;
  orcaref[1] << 0.069946434, 0.000000000, -0.070828201, -0.050226631, -0.000000000, 0.085945300, -0.019711036,
      0.000000000, -0.015109295;
  orcaref[2] << 0.006253782, 0.000000000, -0.121327195, -0.115531271, -0.000000000, 0.037555768, 0.109275799,
      0.000000000, 0.083771206;
  orcaref[3] << 0.065862949, 0.000000000, -0.101464090, -0.080868158, -0.000000000, 0.089958397, 0.015004295,
      0.000000000, 0.011506103;

  for (unsigned i = 0; i < 4; i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-5);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDA_R_LDA_DENSFITJ) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::LDA;
  settings.basis.label = "def2-svp";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.auxJLabel = "def2-svp-ri-c";
  settings.basis.auxCLabel = "def2-svp-ri-c";
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDA;
  lrscf.settings.densFitJ = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 2;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 6.0.1 January 2025
  std::vector<Eigen::MatrixXd> orcaref(2, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = +0.1699300529;
  orcaref[0](1, 2) = -0.1699300529;
  orcaref[1](0, 2) = +0.0318901337;
  orcaref[1](1, 2) = -0.0318901337;

  for (unsigned i = 0; i < orcaref.size(); i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-6);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDDFT_U_GGA_DENSFITJ) {
  Settings settings;
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::PBE;
  settings.basis.label = "def2-svp";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.auxJLabel = "def2-svp-ri-c";
  settings.basis.auxCLabel = "def2-svp-ri-c";
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::UNRESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 2;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2};
  TDDFTGradientCalculator<UNRESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 6.0.1 January 2025
  std::vector<Eigen::MatrixXd> orcaref(2, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = 0.281110667;
  orcaref[0](1, 2) = -0.281110667;
  orcaref[1](0, 2) = 0.178031585;
  orcaref[1](1, 2) = -0.178031585;

  for (unsigned i = 0; i < orcaref.size(); i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 5e-6);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

TEST_F(TDDFTGradientsTest, TDDFT_R_RangeSeparatedHybridGGA_DENSFITJ) {
  Settings settings;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.label = "def2-svp";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.grid.accuracy = 7;
  settings.grid.smallGridAccuracy = 7;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.auxJLabel = "ri_j_weigend";
  settings.basis.auxCLabel = "ri_j_weigend";
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings, 0, 0);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf({system});
  auto lrscfContr = std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(system, lrscf.settings);

  lrscf.settings.method = Options::LR_METHOD::TDDFT;
  lrscf.settings.densFitJ = Options::DENS_FITS::RI;
  lrscf.settings.nEigen = 3;
  lrscf.settings.grid.accuracy = 7;
  lrscf.settings.grid.smallGridAccuracy = 7;
  lrscf.settings.conv = 5e-8;
  lrscf.run();
  lrscf.settings.excGradList = {1, 2, 3};
  TDDFTGradientCalculator<RESTRICTED> tddftgrads(lrscfContr);
  std::vector<Eigen::MatrixXd> excgrads = tddftgrads.calculateGradients();

  // ORCA 6.0.1 January 2025
  std::vector<Eigen::MatrixXd> orcaref(3, Eigen::MatrixXd::Zero(2, 3));
  orcaref[0](0, 2) = +0.180881728;
  orcaref[0](1, 2) = -0.180881728;
  orcaref[1](0, 2) = +0.021617891;
  orcaref[1](1, 2) = -0.021617891;
  orcaref[2](0, 2) = +0.517438909;
  orcaref[2](1, 2) = -0.517438909;

  for (unsigned i = 0; i < orcaref.size(); i++) {
    for (unsigned at = 0; at < orcaref[i].rows(); at++) {
      for (unsigned j = 0; j < 3; j++) {
        EXPECT_NEAR(excgrads[i](at, j), orcaref[i](at, j), 1e-5);
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
}

} // namespace Serenity
