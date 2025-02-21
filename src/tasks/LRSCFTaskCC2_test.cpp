/**
 * @file LRSCFTaskCC2_test.cpp
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

/* Include Serenity Internal Headers */
#include "geometry/AtomTypeFactory.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class LRSCFTaskCC2Test
 * @brief Sets everything up for the tests of RICC2.h/.cpp.\n
          All of these results were obtained with Turbomole 7.4.1 on January 16, 2021.
 */
class LRSCFTaskCC2Test : public ::testing::Test {
 protected:
  LRSCFTaskCC2Test()
    : sys(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ)) {
  }

  virtual ~LRSCFTaskCC2Test() = default;

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
TEST_F(LRSCFTaskCC2Test, CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, SCS_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;

  task.settings.sss = 0.333;
  task.settings.oss = 1.200;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 2e-8;

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
TEST_F(LRSCFTaskCC2Test, NAF_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, SCS_uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, CISDINF) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
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
 * @brief Tests LRSCFTask: CISD
 */
TEST_F(LRSCFTaskCC2Test, CISD) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISD;
  task.settings.cctrdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-6;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  auto compare = [&]() {
    // obtained with turbomole 7.4.1, june 2020

    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.150083089917, accuracy);
    EXPECT_NEAR(results(1, 0), 0.235579424557, accuracy);
    EXPECT_NEAR(results(2, 0), 0.272190287743, accuracy);
    EXPECT_NEAR(results(3, 0), 0.278107905930, accuracy);

    EXPECT_NEAR(results(0, 1), 0.000000000041, accuracy);
    EXPECT_NEAR(results(1, 1), 0.018401805944, accuracy);
    EXPECT_NEAR(results(2, 1), 0.039706677753, accuracy);
    EXPECT_NEAR(results(3, 1), 0.116076563277, accuracy);
  };

  // test with standard settings
  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* CISD */

/**
 * @test
 * @brief Tests LRSCFTask: SCS-CISDINF
 */
TEST_F(LRSCFTaskCC2Test, SCS_CISDINF) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
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
TEST_F(LRSCFTaskCC2Test, uCISDINF) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, SCS_uCISDINF) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CISDINF;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, NAF_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, SCS_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, SCS_uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
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
TEST_F(LRSCFTaskCC2Test, CC2_Chiral) {
  auto diaz = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({diaz});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  Eigen::MatrixXd excitations(task.settings.nEigen, 7);
  excitations.row(0) << 0.229138352906, 0.000436648226, 0.000437489786, -0.000137019112, -0.000137213965,
      -0.000137019112, 0.265233358128;
  excitations.row(1) << 0.238283887674, 0.025289430337, 0.023459798849, -0.000190537608, -0.000179037740,
      -0.000185755524, 0.184595342429;
  excitations.row(2) << 0.256274773647, 0.003381275026, 0.003365279814, -0.000125052299, -0.000124803822,
      -0.000125052299, 1.649061366311;
  excitations.row(3) << 0.261056519323, 0.073754001921, 0.066161789589, 0.000405431775, 0.000371064228, 0.000391520578,
      1.355593205839;

  EXPECT_LE((results - excitations).cwiseAbs().maxCoeff(), accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(diaz);
  SystemController__TEST_SUPPLY::cleanUp();
} /* CC2-Chiral */

/**
 * @test
 * @brief Tests LRSCFTask: ADC2-Chiral
 */
TEST_F(LRSCFTaskCC2Test, ADC2_Chiral) {
  auto diaz = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({diaz});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.nafThresh = 1e-2;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-7;

  Eigen::MatrixXd excitations(task.settings.nEigen, 7);
  excitations.row(0) << 0.227400601818, 0.000313078708, 0.000558341587, -0.000115376764, -0.000154078335,
      -0.000115376764, 0.092393859770;
  excitations.row(1) << 0.236584525153, 0.020715602186, 0.026296601386, -0.000179660127, -0.000195230500,
      -0.000173279280, 0.047133852752;
  excitations.row(2) << 0.254475721277, 0.001969010639, 0.003386148702, -0.000099621707, -0.000130642042,
      -0.000099621707, 1.795190878687;
  excitations.row(3) << 0.259201255789, 0.061986565053, 0.071257291943, 0.000364241819, 0.000365364782, 0.000340769574,
      1.572581491847;

  EXPECT_LE((results - excitations).cwiseAbs().maxCoeff(), accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(diaz);
  SystemController__TEST_SUPPLY::cleanUp();
} /* ADC2-Chiral */

/**
 * @test
 * @brief Tests LRSCFTask: Cholesky (ACCD) ADC2
 */
TEST_F(LRSCFTaskCC2Test, CD_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.densFitCache = Options::DENS_FITS::ACCD;
  task.run();
  auto& results = task.getTransitions();

  double accuracy = 1e-8;

  // Serenity Oct 2021
  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.145392596984, accuracy);
    EXPECT_NEAR(results(1, 0), 0.231148491513, accuracy);
    EXPECT_NEAR(results(2, 0), 0.267311180261, accuracy);
    EXPECT_NEAR(results(3, 0), 0.271244285999, accuracy);

    // oscillator strengths length
    EXPECT_NEAR(results(0, 1), 0.000000000028, accuracy);
    EXPECT_NEAR(results(1, 1), 0.016385152217, accuracy);
    EXPECT_NEAR(results(2, 1), 0.045497875605, accuracy);
    EXPECT_NEAR(results(3, 1), 0.027898292568, accuracy);

    // oscillator strengths velocity
    EXPECT_NEAR(results(0, 2), 0.000000000018, accuracy);
    EXPECT_NEAR(results(1, 2), 0.017752006296, accuracy);
    EXPECT_NEAR(results(2, 2), 0.061060664733, accuracy);
    EXPECT_NEAR(results(3, 2), 0.041156098824, accuracy);
  };

  compare();
} /* ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT-CC2 (approximate embedding, simple FDEc, singlet)
 */
TEST_F(LRSCFTaskCC2Test, FDEc_CC2) {
  // Just use the def2-SVP correlation basis as there is none for Pople bases.
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, true);
  Settings settings_act = act->getSettings();
  Settings settings_env = env->getSettings();
  settings_act.basis.auxCLabel = "DEF2-SVP-RI-C";
  settings_env.basis.auxCLabel = "DEF2-SVP-RI-C";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings_act);
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, settings_env);

  auto task = FreezeAndThawTask<RESTRICTED>({act, env});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.maxCycles = 3;
  task.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.nEigen = 4;
  lrscfA.settings.method = Options::LR_METHOD::CC2;
  lrscfA.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfA.settings.densFitK = Options::DENS_FITS::RI;
  lrscfA.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfA.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  lrscfA.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.nEigen = 4;
  lrscfB.settings.method = Options::LR_METHOD::CC2;
  lrscfB.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfB.settings.densFitK = Options::DENS_FITS::RI;
  lrscfB.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  lrscfB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  lrscfAB.settings.method = Options::LR_METHOD::CC2;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitK = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfAB.settings.ccexdens = true;
  lrscfAB.settings.cctrdens = true;
  lrscfAB.run();

  Eigen::MatrixXd ref(8, 6);
  ref.row(0) << 0.323758316584, 0.017227013642, 0.062371586920, 0.000012762898, 0.000027040678, 0.000014211296;
  ref.row(1) << 0.324950306395, 0.019494108424, 0.068641074537, -0.000014935770, -0.000029185971, -0.000015554206;
  ref.row(2) << 0.404745872689, 0.000240735232, 0.000174081358, 0.000001347055, 0.000001049588, 0.000001234475;
  ref.row(3) << 0.421237967796, 0.073492472018, 0.117937372982, 0.000120630110, 0.000143989184, 0.000113690299;
  ref.row(4) << 0.427786760057, 0.154264670853, 0.222034732416, -0.000126402157, -0.000152660066, -0.000127277700;
  ref.row(5) << 0.430946240290, 0.002281371711, 0.004901382001, 0.000002368816, 0.000002401562, 0.000001638332;
  ref.row(6) << 0.511291734220, 0.096984129949, 0.089449798591, 0.000030181967, 0.000029520825, 0.000030762586;
  ref.row(7) << 0.534716274423, 0.081331902036, 0.074173194419, -0.000021416661, -0.000018764943, -0.000019662657;

  EXPECT_LT((ref.col(0) - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff(), 1e-7);
  EXPECT_LT((ref.col(1) - lrscfAB.getTransitions().col(1)).cwiseAbs().maxCoeff(), 1e-7);
  EXPECT_LT((ref.col(2) - lrscfAB.getTransitions().col(2)).cwiseAbs().maxCoeff(), 1e-7);
  EXPECT_LT((ref.col(3) - lrscfAB.getTransitions().col(3)).cwiseAbs().maxCoeff(), 1e-7);
  EXPECT_LT((ref.col(4) - lrscfAB.getTransitions().col(4)).cwiseAbs().maxCoeff(), 1e-7);
  EXPECT_LT((ref.col(5) - lrscfAB.getTransitions().col(5)).cwiseAbs().maxCoeff(), 1e-7);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LRSCFTask: ADC2 (exact embedding, full FDEc, triplet)
 */
TEST_F(LRSCFTaskCC2Test, FDEc_ADC2) {
  // Just use the def2-SVP correlation basis as there is none for Pople bases.
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, true);
  Settings settings_act = act->getSettings();
  Settings settings_env = env->getSettings();
  settings_act.basis.auxCLabel = "DEF2-SVP-RI-C";
  settings_env.basis.auxCLabel = "DEF2-SVP-RI-C";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings_act);
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, settings_env);

  auto task = TDEmbeddingTask<RESTRICTED>(act, env);
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  task.run();

  auto task2 = FDETask<RESTRICTED>(env, {act});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  task2.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {env});
  lrscfA.settings.nEigen = 4;
  lrscfA.settings.method = Options::LR_METHOD::ADC2;
  lrscfA.settings.excludeProjection = true;
  lrscfA.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfA.settings.densFitK = Options::DENS_FITS::RI;
  lrscfA.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfA.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  lrscfA.settings.triplet = true;
  lrscfA.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfB({env}, {act});
  lrscfB.settings.nEigen = 4;
  lrscfB.settings.method = Options::LR_METHOD::ADC2;
  lrscfB.settings.excludeProjection = true;
  lrscfB.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfB.settings.densFitK = Options::DENS_FITS::RI;
  lrscfB.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  lrscfB.settings.triplet = true;
  lrscfB.run();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  lrscfAB.settings.method = Options::LR_METHOD::ADC2;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitK = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::HF;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  lrscfAB.settings.triplet = true;
  lrscfAB.run();

  Eigen::MatrixXd ref(8, 1);
  ref.col(0) << 0.296354025192, 0.296574204636, 0.379028220743, 0.386474296341, 0.388211587389, 0.393771155946,
      0.462816909506, 0.474120969096;
  EXPECT_LT((ref.col(0) - lrscfAB.getTransitions().col(0)).cwiseAbs().maxCoeff(), 1e-7);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LRSCFTask: LT-CC2
 */
#ifdef SERENITY_USE_LAPLACE_MINIMAX
TEST_F(LRSCFTaskCC2Test, LT_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::CC2;
  tasklt.settings.cctrdens = true;
  tasklt.settings.ccexdens = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-12;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  double accuracy = 1e-7;

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT_uCC2
 */
TEST_F(LRSCFTaskCC2Test, LT_uCC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);

  double accuracy = 1e-7;

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::CC2;
  tasklt.settings.cctrdens = true;
  tasklt.settings.ccexdens = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-12;
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
TEST_F(LRSCFTaskCC2Test, LT_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::RESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::ADC2;
  tasklt.settings.cctrdens = true;
  tasklt.settings.ccexdens = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-12;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  double accuracy = 1e-7;

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-ADC2 */

/**
 * @test
 * @brief Tests LRSCFTask: uADC2
 */
TEST_F(LRSCFTaskCC2Test, LT_uADC2) {
  sys->setSpin(1);
  sys->setCharge(1);
  sys->setSCFMode(Options::SCF_MODES::UNRESTRICTED);

  double accuracy = 1e-7;

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.sss = 0.0;
  task.settings.oss = 1.3;
  task.run();
  auto& results = task.getTransitions();

  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> tasklt({sys});
  tasklt.settings.method = Options::LR_METHOD::ADC2;
  tasklt.settings.cctrdens = true;
  tasklt.settings.ccexdens = true;
  tasklt.settings.preopt = 1e-6;
  tasklt.settings.conv = 1e-9;
  tasklt.settings.ltconv = 1e-13;
  tasklt.settings.sss = 0.0;
  tasklt.settings.oss = 1.3;
  tasklt.run();
  auto& resultslt = tasklt.getTransitions();

  double maxDiff = (results - resultslt).cwiseAbs().maxCoeff();

  EXPECT_LT(maxDiff, accuracy);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* LT-uCC2 */

TEST_F(LRSCFTaskCC2Test, Response_Properties_CC2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ);
  auto settings = sys->getSettings();

  // Strict convergence criteria.
  settings.scf.energyThreshold = 1e-9;
  settings.scf.diisThreshold = 1e-9;
  settings.scf.rmsdThreshold = 1e-9;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Diaziridine_HF_AUG_CC_PVDZ, settings);

  // Run in undamped mode.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf1({sys});
  lrscf1.settings.nEigen = 4;
  lrscf1.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf1.settings.method = Options::LR_METHOD::CC2;
  lrscf1.settings.cctrdens = true;
  lrscf1.settings.ccexdens = true;
  lrscf1.settings.frozenCore = true;
  lrscf1.settings.nafThresh = 1e-2;
  lrscf1.settings.conv = 1e-7;
  lrscf1.settings.diis = true;
  lrscf1.run();

  // Run in undamped mode (unrestricted).
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf2({sys});
  lrscf2.settings.nEigen = 4;
  lrscf2.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf2.settings.method = Options::LR_METHOD::CC2;
  lrscf2.settings.cctrdens = true;
  lrscf2.settings.ccexdens = true;
  lrscf2.settings.frozenCore = true;
  lrscf2.settings.nafThresh = 1e-2;
  lrscf2.settings.conv = 1e-7;
  lrscf2.settings.diis = true;
  lrscf2.run();

  // Run in undamped mode (velocity gauge).
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf3({sys});
  lrscf3.settings.nEigen = 4;
  lrscf3.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf3.settings.gauge = Options::GAUGE::VELOCITY;
  lrscf3.settings.method = Options::LR_METHOD::CC2;
  lrscf3.settings.cctrdens = true;
  lrscf3.settings.ccexdens = true;
  lrscf3.settings.frozenCore = true;
  lrscf3.settings.nafThresh = 1e-2;
  lrscf3.settings.conv = 1e-7;
  lrscf3.settings.diis = true;
  lrscf3.run();

  // Run in undamped mode (velocity gauge, unrestricted).
  LRSCFTask<Options::SCF_MODES::UNRESTRICTED> lrscf4({sys});
  lrscf4.settings.nEigen = 4;
  lrscf4.settings.frequencies = {3.0, 4.0, 5.0};
  lrscf4.settings.gauge = Options::GAUGE::VELOCITY;
  lrscf4.settings.method = Options::LR_METHOD::CC2;
  lrscf4.settings.cctrdens = true;
  lrscf4.settings.ccexdens = true;
  lrscf4.settings.frozenCore = true;
  lrscf4.settings.nafThresh = 1e-2;
  lrscf4.settings.conv = 1e-7;
  lrscf4.settings.diis = true;
  lrscf4.run();

  Eigen::MatrixXd excitations_restricted(lrscf1.settings.nEigen, 7);
  Eigen::MatrixXd excitations_unrestricted(lrscf1.settings.nEigen, 7);
  excitations_restricted.row(0) << 0.229154658748, 0.000434129156, 0.000426134684, -0.000137002911, -0.000135800501,
      -0.000137002911, 0.294810214007;
  excitations_restricted.row(1) << 0.238278316311, 0.025347213746, 0.024395468471, -0.000190445262, -0.000181559423,
      -0.000184936485, 0.213763889860;
  excitations_restricted.row(2) << 0.256300777641, 0.003398755066, 0.003403582643, -0.000126085965, -0.000126227205,
      -0.000126085965, 1.616200611912;
  excitations_restricted.row(3) << 0.261083541956, 0.073825834241, 0.066876352482, 0.000408479366, 0.000374976073,
      0.000393725257, 1.321923336540;

  excitations_unrestricted.row(0) << 0.224687509373, 0.000000000000, 0.000000000000, -0.000000000000, -0.000000000000,
      -0.000000000000, 0.580727850969;
  excitations_unrestricted.row(1) << 0.229154654354, 0.000434134318, 0.000426139068, -0.000137003852, -0.000135801389,
      -0.000137003852, 0.294815594623;
  excitations_unrestricted.row(2) << 0.231154473407, 0.000000000000, 0.000000000001, -0.000000000000, -0.000000000000,
      -0.000000000000, 0.266076668725;
  excitations_unrestricted.row(3) << 0.255847406924, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
      0.000000000000, 0.900372242445;

  Eigen::MatrixXd undamped_properties(lrscf1.settings.frequencies.size(), 5);
  undamped_properties.row(0) << 0.110247967, 32.5925094013, -0.1765344140, 0.000000000, 0.000000000;
  undamped_properties.row(1) << 0.146997290, 34.2005005175, -0.2746899944, 0.000000000, 0.000000000;
  undamped_properties.row(2) << 0.183746612, 36.9506338012, -0.5926584903, 0.000000000, 0.000000000;

  Eigen::MatrixXd undamped_properties_velocity(lrscf1.settings.frequencies.size(), 5);
  undamped_properties_velocity.row(0) << 0.110247967, 29.9122759486, -0.2141005946, 0.000000000, 0.000000000;
  undamped_properties_velocity.row(1) << 0.146997290, 31.3894495407, -0.3194350413, 0.000000000, 0.000000000;
  undamped_properties_velocity.row(2) << 0.183746612, 33.9213362978, -0.6500083749, 0.000000000, 0.000000000;

  // Compare excitation energies and oscillator/rotator strengths.
  EXPECT_LE((excitations_restricted - lrscf1.getTransitions()).cwiseAbs().maxCoeff(), 1e-6);
  std::cout << "unrestricted reference is\n" << excitations_unrestricted << std::endl;
  std::cout << "unres result in length gauge is\n" << lrscf2.getTransitions() << std::endl;
  std::cout << "their difference is\n" << (excitations_unrestricted - lrscf4.getTransitions()).cwiseAbs() << std::endl;
  EXPECT_LE((excitations_unrestricted - lrscf2.getTransitions()).cwiseAbs().maxCoeff(), 1e-6);
  EXPECT_LE((excitations_restricted - lrscf3.getTransitions()).cwiseAbs().maxCoeff(), 1e-6);
  std::cout << "unres result in velocity gauge is\n" << lrscf4.getTransitions() << std::endl;
  std::cout << "their difference is\n" << (excitations_unrestricted - lrscf4.getTransitions()).cwiseAbs() << std::endl;
  EXPECT_LE((excitations_unrestricted - lrscf4.getTransitions()).cwiseAbs().maxCoeff(), 1e-6);

  // Compare undamped properties.
  for (unsigned iFreq = 0; iFreq < 3; ++iFreq) {
    double frequency = lrscf1.settings.frequencies[iFreq];
    double isoPol1 = 1. / 3. * std::get<1>(lrscf1.getProperties()[iFreq]).trace();
    double isoPol2 = 1. / 3. * std::get<1>(lrscf2.getProperties()[iFreq]).trace();
    double isoOR1 = -1. / 3. / frequency * std::get<2>(lrscf1.getProperties()[iFreq]).trace();
    double isoOR2 = -1. / 3. / frequency * std::get<2>(lrscf2.getProperties()[iFreq]).trace();
    double isoAbs1 = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf1.getProperties()[iFreq]).trace();
    double isoAbs2 = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf2.getProperties()[iFreq]).trace();
    double isoCD1 = 6.533 * frequency * std::get<4>(lrscf1.getProperties()[iFreq]).trace();
    double isoCD2 = 6.533 * frequency * std::get<4>(lrscf2.getProperties()[iFreq]).trace();
    EXPECT_NEAR(isoPol1, undamped_properties(iFreq, 1), 1E-6);
    EXPECT_NEAR(isoPol2, undamped_properties(iFreq, 1), 1E-6);
    EXPECT_NEAR(isoOR1, undamped_properties(iFreq, 2), 1E-6);
    EXPECT_NEAR(isoOR2, undamped_properties(iFreq, 2), 1E-6);
    EXPECT_NEAR(isoAbs1, undamped_properties(iFreq, 3), 1E-6);
    EXPECT_NEAR(isoAbs2, undamped_properties(iFreq, 3), 1E-6);
    EXPECT_NEAR(isoCD1, undamped_properties(iFreq, 4), 1E-6);
    EXPECT_NEAR(isoCD2, undamped_properties(iFreq, 4), 1E-6);
  }
  // Compare undamped properties (velocity).
  for (unsigned iFreq = 0; iFreq < 3; ++iFreq) {
    double frequency = lrscf3.settings.frequencies[iFreq];
    double isoPol1 = 1. / 3. * std::get<1>(lrscf3.getProperties()[iFreq]).trace();
    double isoPol2 = 1. / 3. * std::get<1>(lrscf4.getProperties()[iFreq]).trace();
    double isoOR1 = -1. / 3. / frequency * std::get<2>(lrscf3.getProperties()[iFreq]).trace();
    double isoOR2 = -1. / 3. / frequency * std::get<2>(lrscf4.getProperties()[iFreq]).trace();
    double isoAbs1 = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf3.getProperties()[iFreq]).trace();
    double isoAbs2 = 4. * PI * frequency / 3. / SPEEDOFLIGHT_AU * std::get<3>(lrscf4.getProperties()[iFreq]).trace();
    double isoCD1 = 6.533 * frequency * std::get<4>(lrscf3.getProperties()[iFreq]).trace();
    double isoCD2 = 6.533 * frequency * std::get<4>(lrscf4.getProperties()[iFreq]).trace();
    EXPECT_NEAR(isoPol1, undamped_properties_velocity(iFreq, 1), 1E-6);
    EXPECT_NEAR(isoPol2, undamped_properties_velocity(iFreq, 1), 1E-6);
    EXPECT_NEAR(isoOR1, undamped_properties_velocity(iFreq, 2), 1E-6);
    EXPECT_NEAR(isoOR2, undamped_properties_velocity(iFreq, 2), 1E-6);
    EXPECT_NEAR(isoAbs1, undamped_properties_velocity(iFreq, 3), 1E-6);
    EXPECT_NEAR(isoAbs2, undamped_properties_velocity(iFreq, 3), 1E-6);
    EXPECT_NEAR(isoCD1, undamped_properties_velocity(iFreq, 4), 1E-6);
    EXPECT_NEAR(isoCD2, undamped_properties_velocity(iFreq, 4), 1E-6);
  }

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* Response Properties CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: LT-CC2 (approximate embedding, response properties)
 */
TEST_F(LRSCFTaskCC2Test, FDEc_CC2_Prop) {
  // Just use the def2-SVP correlation basis as there is none for Pople bases.
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, true);
  Settings settings_act = act->getSettings();
  Settings settings_env = env->getSettings();
  settings_act.basis.auxCLabel = "DEF2-SVP-RI-C";
  settings_env.basis.auxCLabel = "DEF2-SVP-RI-C";
  act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings_act);
  env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs, settings_env);

  auto task = FreezeAndThawTask<RESTRICTED>({act, env});
  task.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task.settings.maxCycles = 3;
  task.run();

  // Actually coupled calculation.
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAB({act, env}, {});
  lrscfAB.settings.method = Options::LR_METHOD::CC2;
  lrscfAB.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitK = Options::DENS_FITS::RI;
  lrscfAB.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfAB.settings.frequencies = {0};
  lrscfAB.settings.nEigen = 0;
  lrscfAB.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  lrscfAB.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  lrscfAB.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  lrscfAB.run();
  EXPECT_NEAR(1.0 / 3.0 * std::get<1>(lrscfAB.getProperties()[0]).trace(), 9.970283807714, 1e-8);

  // Sanity Check: Two times the same subsystem without any coupling
  // should yield two times the isolated polarizability
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfA({act}, {});
  lrscfA.settings.method = Options::LR_METHOD::CC2;
  lrscfA.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfA.settings.densFitK = Options::DENS_FITS::RI;
  lrscfA.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfA.settings.frequencies = {0};
  lrscfA.settings.nEigen = 0;
  lrscfA.run();
  EXPECT_NEAR(1.0 / 3.0 * std::get<1>(lrscfA.getProperties()[0]).trace(), 4.976402942018, 1e-8);

  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscfAA({act, act}, {});
  lrscfAA.settings.method = Options::LR_METHOD::CC2;
  lrscfAA.settings.densFitJ = Options::DENS_FITS::RI;
  lrscfAA.settings.densFitK = Options::DENS_FITS::RI;
  lrscfAA.settings.densFitCache = Options::DENS_FITS::RI;
  lrscfAA.settings.noCoupling = true;
  lrscfAA.settings.frequencies = {0};
  lrscfAA.settings.nEigen = 0;
  lrscfAA.run();
  EXPECT_NEAR(2 * std::get<1>(lrscfA.getProperties()[0]).trace(), std::get<1>(lrscfAA.getProperties()[0]).trace(), 1e-8);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(act);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests LRSCFTask: Triplet-SOS-CC2
 */
TEST_F(LRSCFTaskCC2Test, Triplet_SOS_CC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.triplet = true;

  task.settings.sss = 0.000;
  task.settings.oss = 1.300;
  task.run();
  auto results = task.getTransitions();

  double accuracy = 1e-6;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1410557, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2294341, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2518746, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2871732, accuracy);
  };

  compare();

  task = LRSCFTask<Options::SCF_MODES::RESTRICTED>({sys});
  task.settings.method = Options::LR_METHOD::CC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.triplet = true;
  task.settings.ltconv = 1e-8;

  task.settings.sss = 0.000;
  task.settings.oss = 1.300;
  task.run();
  results = task.getTransitions();

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* Triplet-SOS-CC2 */

/**
 * @test
 * @brief Tests LRSCFTask: Triplet-SOS-ADC2
 */
TEST_F(LRSCFTaskCC2Test, Triplet_SOS_ADC2) {
  LRSCFTask<Options::SCF_MODES::RESTRICTED> task({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.triplet = true;

  task.settings.sss = 0.000;
  task.settings.oss = 1.300;
  task.run();
  auto results = task.getTransitions();

  double accuracy = 1e-6;

  auto compare = [&]() {
    // excitation energies
    EXPECT_NEAR(results(0, 0), 0.1361518, accuracy);
    EXPECT_NEAR(results(1, 0), 0.2266177, accuracy);
    EXPECT_NEAR(results(2, 0), 0.2476206, accuracy);
    EXPECT_NEAR(results(3, 0), 0.2833782, accuracy);
  };

  compare();

  task = LRSCFTask<Options::SCF_MODES::RESTRICTED>({sys});
  task.settings.method = Options::LR_METHOD::ADC2;
  task.settings.cctrdens = true;
  task.settings.ccexdens = true;
  task.settings.preopt = 1e-6;
  task.settings.conv = 1e-9;
  task.settings.triplet = true;
  task.settings.ltconv = 1e-8;

  task.settings.sss = 0.000;
  task.settings.oss = 1.300;
  task.run();
  results = task.getTransitions();

  compare();

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
  SystemController__TEST_SUPPLY::cleanUp();
} /* Triplet-SOS-ADC2 */
#endif /* SERENITY_USE_LAPLACE_MINIMAX */

} // namespace Serenity
