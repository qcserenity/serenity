/**
 * @file LRSCFTaskXWF_test.cpp
 *
 *  @date      May 29, 2020
 *  @author    Niklas Niemeyer
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

#include "geometry/AtomTypeFactory.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFTaskXWFTest
 * @brief Sets everything up for the tests of RICC2.h/.cpp.\n
          All of these results were obtained with Turbomole 7.4.1 on January 16, 2021.
 */
class LRSCFTaskXWFTest : public ::testing::Test {
 protected:
  LRSCFTaskXWFTest()
    : sys(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ)) {
  }

  virtual ~LRSCFTaskXWFTest() = default;

  /// system
  std::shared_ptr<SystemController> sys;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests LRSCFTask: CC2
 */
TEST_F(LRSCFTaskXWFTest, CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1515258005, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2355965653, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2716777268, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2751246394, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.02125456, accuracy);
    EXPECT_NEAR(results(2, 1), 0.05907212, accuracy);
    EXPECT_NEAR(results(3, 1), 0.03609525, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.01846475, accuracy);
    EXPECT_NEAR(results(2, 2), 0.05663694, accuracy);
    EXPECT_NEAR(results(3, 2), 0.03555224, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-CC2
 */
TEST_F(LRSCFTaskXWFTest, SCS_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1552202543, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2487225084, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2844927080, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2861883668, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.01906270, accuracy);
    EXPECT_NEAR(results(2, 1), 0.05958140, accuracy);
    EXPECT_NEAR(results(3, 1), 0.04087179, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.01718568, accuracy);
    EXPECT_NEAR(results(2, 2), 0.05652130, accuracy);
    EXPECT_NEAR(results(3, 2), 0.04169574, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: NAF-CC2
 */
TEST_F(LRSCFTaskXWFTest, NAF_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.151349862, accuracy);
    EXPECT_NEAR(results(1, 0), 0.235674588, accuracy);
    EXPECT_NEAR(results(2, 0), 0.271868330, accuracy);
    EXPECT_NEAR(results(3, 0), 0.275291449, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.021415750, accuracy);
    EXPECT_NEAR(results(2, 1), 0.059244388, accuracy);
    EXPECT_NEAR(results(3, 1), 0.035790577, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.018535595, accuracy);
    EXPECT_NEAR(results(2, 2), 0.056819540, accuracy);
    EXPECT_NEAR(results(3, 2), 0.035396143, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* NAF-CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: uCC2
 */
TEST_F(LRSCFTaskXWFTest, uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1465144776, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2017040472, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2281529165, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2479700517, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.00031056, accuracy);
    EXPECT_NEAR(results(2, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 1), 0.02673360, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.00066733, accuracy);
    EXPECT_NEAR(results(2, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 2), 0.01829425, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* uCC2 */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-uCC2
 */
TEST_F(LRSCFTaskXWFTest, SCS_uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1421235742, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2004130417, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2307112080, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2552132338, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.00057084, accuracy);
    EXPECT_NEAR(results(2, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 1), 0.02774498, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.00124519, accuracy);
    EXPECT_NEAR(results(2, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 2), 0.01610219, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-uCC2 */

/**
 * @test
 * @brief Tests LRSCFTask: CISDINF
 */
TEST_F(LRSCFTaskXWFTest, CISDINF) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.ccprops = true;
  task.settings.conv = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // obtained with turbomole 7.4.1, june 2020

    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1453855, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2311446, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2673037, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2712229, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* CISDINF */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-CISDINF
 */
TEST_F(LRSCFTaskXWFTest, SCS_CISDINF) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.ccprops = true;
  task.settings.conv = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // obtained with turbomole 7.4.1, june 2020

    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1494589, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2443141, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2801659, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2823102, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-CISDINF */

/**
 * @test
 * @brief Tests LRSCFTask: uCISDINF
 */
TEST_F(LRSCFTaskXWFTest, uCISDINF) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1418819, accuracy);
    EXPECT_NEAR(results(1, 0), 0.1983554, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2256536, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2581843, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* uCISDINF */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-uCISDINF
 */
TEST_F(LRSCFTaskXWFTest, SCS_uCISDINF) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1395063, accuracy);
    EXPECT_NEAR(results(1, 0), 0.1984721, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2287448, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2703374, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-uCISDINF */

/**
 * @test
 * @brief Tests LRSCFTask: ADC2
 */
TEST_F(LRSCFTaskXWFTest, ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1453638619, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2311434937, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2673026271, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2712218925, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.01639469, accuracy);
    EXPECT_NEAR(results(2, 1), 0.04550647, accuracy);
    EXPECT_NEAR(results(3, 1), 0.02789352, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.01776513, accuracy);
    EXPECT_NEAR(results(2, 2), 0.06108359, accuracy);
    EXPECT_NEAR(results(3, 2), 0.04115956, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: NAF-ADC2
 */
TEST_F(LRSCFTaskXWFTest, NAF_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.145185360, accuracy);
    EXPECT_NEAR(results(1, 0), 0.231198642, accuracy);
    EXPECT_NEAR(results(2, 0), 0.267465763, accuracy);
    EXPECT_NEAR(results(3, 0), 0.271374929, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.016518787, accuracy);
    EXPECT_NEAR(results(2, 1), 0.045662282, accuracy);
    EXPECT_NEAR(results(3, 1), 0.027658924, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.017833051, accuracy);
    EXPECT_NEAR(results(2, 2), 0.061295481, accuracy);
    EXPECT_NEAR(results(3, 2), 0.041043368, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* NAF-ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-ADC2
 */
TEST_F(LRSCFTaskXWFTest, SCS_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1494353469, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2443130213, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2801647240, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2823091303, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.01599554, accuracy);
    EXPECT_NEAR(results(2, 1), 0.04840635, accuracy);
    EXPECT_NEAR(results(3, 1), 0.03239626, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.01658132, accuracy);
    EXPECT_NEAR(results(2, 2), 0.05992600, accuracy);
    EXPECT_NEAR(results(3, 2), 0.04800029, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: uADC2
 */
TEST_F(LRSCFTaskXWFTest, uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1418804232, accuracy);
    EXPECT_NEAR(results(1, 0), 0.1983517132, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2256431422, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2581787577, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.00123054, accuracy);
    EXPECT_NEAR(results(2, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 1), 0.01486206, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.00073673, accuracy);
    EXPECT_NEAR(results(2, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 2), 0.02182078, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* uADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-uADC2
 */
TEST_F(LRSCFTaskXWFTest, SCS_uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1395039154, accuracy);
    EXPECT_NEAR(results(1, 0), 0.1984681276, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2287346301, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2702055609, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.00140322, accuracy);
    EXPECT_NEAR(results(2, 1), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 1), 0.02033615, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.00144423, accuracy);
    EXPECT_NEAR(results(2, 2), 0.00000000, accuracy);
    EXPECT_NEAR(results(3, 2), 0.01336291, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* SCS-uADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: CC2-Chiral
 */
TEST_F(LRSCFTaskXWFTest, CC2_Chiral) {
  auto diaz = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({diaz});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  Eigen::MatrixXd excitations(task.settings.nEigen, 6);
  excitations.row(0) << 0.229138353, 0.000436648, 0.000437490, -0.000137019, -0.000137214, -0.000137019;
  excitations.row(1) << 0.238283888, 0.025289430, 0.023459799, -0.000190538, -0.000179038, -0.000185756;
  excitations.row(2) << 0.256274774, 0.003381275, 0.003365280, -0.000125052, -0.000124804, -0.000125052;
  excitations.row(3) << 0.261056519, 0.073754001, 0.066161789, 0.000405432, 0.000371064, 0.000391521;

  EXPECT_LE((results - excitations).cwiseAbs().maxCoeff(), accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(diaz);
  SystemController__TEST_SUPPLY::cleanUp();
} /* CC2-Chiral */

/**
 * @test
 * @brief Tests LRSCFTask: ADC2-Chiral
 */
TEST_F(LRSCFTaskXWFTest, ADC2_Chiral) {
  auto diaz = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({diaz});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  Eigen::MatrixXd excitations(task.settings.nEigen, 6);
  excitations.row(0) << 0.227400602, 0.000313079, 0.000558341, -0.000115377, -0.000154078, -0.000115377;
  excitations.row(1) << 0.236584525, 0.020715602, 0.026296601, -0.000179660, -0.000195230, -0.000173279;
  excitations.row(2) << 0.254475721, 0.001969011, 0.003386149, -0.000099622, -0.000130642, -0.000099622;
  excitations.row(3) << 0.259201256, 0.061986566, 0.071257293, 0.000364242, 0.000365365, 0.000340770;

  EXPECT_LE((results - excitations).cwiseAbs().maxCoeff(), accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(diaz);
  SystemController__TEST_SUPPLY::cleanUp();
} /* ADC2-Chiral */

/**
 * @test
 * @brief Tests LRSCFTask: Cholesky (ACCD) ADC2
 */
TEST_F(LRSCFTaskXWFTest, CD_ADC2) {
  Settings settings = sys->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::ACCD;
  settings.basis.densFitK = Options::DENS_FITS::ACCD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, settings);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.densFitCache = Options::DENS_FITS::ACCD;
  task.settings.densFitJ = Options::DENS_FITS::ACCD;
  task.settings.densFitK = Options::DENS_FITS::ACCD;
  task.settings.densFitLRK = Options::DENS_FITS::ACCD;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  // Serenity Oct 2021
  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.145398597, accuracy);
    EXPECT_NEAR(results(1, 0), 0.231151033, accuracy);
    EXPECT_NEAR(results(2, 0), 0.267305982, accuracy);
    EXPECT_NEAR(results(3, 0), 0.271246114, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 1), 0.016382552, accuracy);
    EXPECT_NEAR(results(2, 1), 0.045487757, accuracy);
    EXPECT_NEAR(results(3, 1), 0.027905454, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.000000000, accuracy);
    EXPECT_NEAR(results(1, 2), 0.017748262, accuracy);
    EXPECT_NEAR(results(2, 2), 0.061045723, accuracy);
    EXPECT_NEAR(results(3, 2), 0.041160018, accuracy);
  };

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT-CC2
 */
#ifdef SERENITY_USE_LAPLACE_MINIMAX
TEST_F(LRSCFTaskXWFTest, LT_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::CC2;
  tasklt.settings.ccprops = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-10;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  double accuracy = 1e-8;

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT_uCC2
 */
TEST_F(LRSCFTaskXWFTest, LT_uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);

  double accuracy = 1e-7;

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::CC2;
  tasklt.settings.ccprops = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-10;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-uCC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT-ADC2
 */
TEST_F(LRSCFTaskXWFTest, LT_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::ADC2;
  tasklt.settings.ccprops = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-10;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  double accuracy = 1e-8;

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: uADC2
 */
TEST_F(LRSCFTaskXWFTest, LT_uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);

  double accuracy = 1e-7;

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.ccprops = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::ADC2;
  tasklt.settings.ccprops = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-10;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-uCC2 */
#endif /* SERENITY_USE_LAPLACE_MINIMAX */

} // namespace Serenity
