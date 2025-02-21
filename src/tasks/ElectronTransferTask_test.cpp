/**
 * @file ElectronTransferTask_test.cpp
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
#include "tasks/ElectronTransferTask.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "postHF/ET/FDEDiabController.h"
#include "postHF/ET/FDEETCalculator.h"
#include "postHF/ET/FDEETController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class ElectronTransferTaskTest
 * @brief Sets everything up for the tests of ElectronTransferTask.h/.cpp .
 */
class ElectronTransferTaskTest : public ::testing::Test {
 protected:
  ElectronTransferTaskTest() {
  }

  virtual ~ElectronTransferTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests joint functionality of ElectronTransferTask.h/.cpp.
 */
TEST_F(ElectronTransferTaskTest, unrestrictedJoint) {
  // get Systems
  auto sys1a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  auto sys2a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  auto sys1b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  auto sys2b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // perform FeezeAndThawTask
  // state1
  auto task1 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1a, sys2a}, {});
  task1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task1.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task1.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task1.settings.maxCycles = 3;
  task1.run();

  // state2
  auto task2 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1b, sys2b}, {});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task2.settings.maxCycles = 3;
  task2.run();

  // perform ElectronTransferTask
  auto task = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  task.settings.states = {2, 2};
  task.settings.couple = {{1, 2}};
  task.settings.spindensity = true;
  task.settings.spinpopulation = true;
  task.settings.population = {0, 1};
  task.run();

  // get values calculated in task
  auto calculator = task.getFDEETCalculator();
  auto H = calculator->getHamiltonian();
  auto S = calculator->getDeterminants();
  auto C = calculator->getLinCoeffs();
  auto E = calculator->getEigenValues();

  // check FDE-ET values
  EXPECT_NEAR(0.2205875, calculator->getAnalyticalCoupling(), 1e-5); // delta E
  EXPECT_NEAR(0.5058070, calculator->getLongRangeExEnergy(), 1e-5);  // V
  EXPECT_NEAR(-15.46291, H(0, 0), 1e-4);                             // H11
  EXPECT_NEAR(-15.29135, H(1, 1), 1e-4);                             // H22
  EXPECT_NEAR(-15.52435, H(0, 1), 1e-4);                             // H12
  EXPECT_NEAR(0.90410, S(0, 0), 1e-4);                               // S11
  EXPECT_NEAR(0.88999, S(1, 1), 1e-4);                               // S22
  EXPECT_NEAR(0.64632, S(0, 1), 1e-4);                               // S12
  EXPECT_NEAR(-0.90515, C(0, 0), 1e-4);                              // C11
  EXPECT_NEAR(1.22585, C(0, 1), 1e-4);                               // C12
  EXPECT_NEAR(-0.19320, C(1, 0), 1e-4);                              // C21
  EXPECT_NEAR(-1.51437, C(1, 1), 1e-4);                              // C22
  EXPECT_NEAR(-15.47110, E(0), 1e-4);                                // E1
  EXPECT_NEAR(-14.96529, E(1), 1e-4);                                // E2

  auto diabController = task.getFDEDiabController();
  auto spinPopsState = diabController->getSortedAtomPopulations();
  // check spin populations
  EXPECT_NEAR(0.69133, spinPopsState[0][0], 1e-4);
  EXPECT_NEAR(0.15424, spinPopsState[0][1], 1e-4);
  EXPECT_NEAR(0.15424, spinPopsState[0][2], 1e-4);
  EXPECT_NEAR(0.44077, spinPopsState[1][0], 1e-4);
  EXPECT_NEAR(0.27962488, spinPopsState[1][1], 1e-4);
  EXPECT_NEAR(0.27962488, spinPopsState[1][2], 1e-4);

  // check <S*S> values
  EXPECT_NEAR(0.78494, diabController->getS2(0), 1e-4);
  EXPECT_NEAR(1.00187, diabController->getS2(1), 1e-4);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // clean up cubes
  std::string basePath = sys1a->getSystemPath().substr(0, sys1a->getSystemPath().size() - sys1a->getSystemName().size() - 1);
  for (unsigned iState = 0; iState < 2; ++iState) {
    std::remove((basePath + "sys-superSystem12/adiabState" + std::to_string(iState) + "_SpinDensity.cube").c_str());
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(basePath + "sys-superSystem12/", "superSystem12");
}
/**
 * @test
 * @brief Tests use of ALMO with GGA functional in ElectronTransferTask.h/.cpp.
 */
TEST_F(ElectronTransferTaskTest, almo_gga) {
  // get Systems
  auto sys1a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  auto sys2a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  auto sys1b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  auto sys2b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // perform FeezeAndThawTask
  // state1
  auto task1 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1a, sys2a});
  task1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::ALMO;
  task1.settings.maxCycles = 3;
  task1.run();

  // state2
  auto task2 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1b, sys2b});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::ALMO;
  task2.settings.maxCycles = 3;
  task2.run();

  // perform ElectronTransferTask
  auto task = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  task.settings.states = {2, 2};
  task.settings.couple = {{1, 2}};
  task.run();

  // get values calculated in task
  auto calculator = task.getFDEETCalculator();
  auto H = calculator->getHamiltonian();
  auto S = calculator->getDeterminants();
  auto C = calculator->getLinCoeffs();
  auto E = calculator->getEigenValues();

  // check FDE-ET values
  EXPECT_NEAR(0.229901, calculator->getAnalyticalCoupling(), 1e-5); // delta E
  EXPECT_NEAR(0.520690, calculator->getLongRangeExEnergy(), 1e-5);  // V
  EXPECT_NEAR(-15.46668, H(0, 0), 1e-4);                            // H11
  EXPECT_NEAR(-15.29857, H(1, 1), 1e-4);                            // H22
  EXPECT_NEAR(-15.53260, H(0, 1), 1e-4);                            // H12
  EXPECT_NEAR(1.05323, S(0, 0), 1e-4);                              // S11
  EXPECT_NEAR(1.04815, S(1, 1), 1e-4);                              // S22
  EXPECT_NEAR(0.76246, S(0, 1), 1e-4);                              // S12
  EXPECT_NEAR(-0.82734, C(0, 0), 1e-4);                             // C11
  EXPECT_NEAR(1.15197, C(0, 1), 1e-4);                              // C12
  EXPECT_NEAR(-0.19092, C(1, 0), 1e-4);                             // C21
  EXPECT_NEAR(-1.40614, C(1, 1), 1e-4);                             // C22
  EXPECT_NEAR(-15.47614, E(0), 1e-4);                               // E1
  EXPECT_NEAR(-14.95545, E(1), 1e-4);                               // E2

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys1a->getSystemPath() + sys1a->getSystemName() + "_tmp/",
                                                        sys1a->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys2a->getSystemPath() + sys2a->getSystemName() + "_tmp/",
                                                        sys2a->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys1b->getSystemPath() + sys1b->getSystemName() + "_tmp/",
                                                        sys1b->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys2b->getSystemPath() + sys2b->getSystemName() + "_tmp/",
                                                        sys2b->getSystemName() + "_tmp");
  std::string basePath = sys1a->getSystemPath().substr(0, sys1a->getSystemPath().size() - sys1a->getSystemName().size() - 1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(basePath + "sys-superSystem12/", "superSystem12");
}
/**
 * @test
 * @brief Tests use of ALMO with hybrid XC functional in ElectronTransferTask.h/.cpp.
 */
TEST_F(ElectronTransferTaskTest, almo_hybrid) {
  // get Systems
  auto sys1a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  auto sys2a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  auto sys1b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  auto sys2b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);
  sys1a->setXCfunctional(CompositeFunctionals::XCFUNCTIONALS::BPW91);
  sys2a->setXCfunctional(CompositeFunctionals::XCFUNCTIONALS::BPW91);
  sys1b->setXCfunctional(CompositeFunctionals::XCFUNCTIONALS::BPW91);
  sys2b->setXCfunctional(CompositeFunctionals::XCFUNCTIONALS::BPW91);

  // perform FeezeAndThawTask
  // state1
  auto task1 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1a, sys2a}, {});
  task1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::ALMO;
  task1.settings.maxCycles = 3;
  task1.run();

  // state2
  auto task2 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1b, sys2b}, {});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::ALMO;
  task2.settings.maxCycles = 3;
  task2.run();

  // perform ElectronTransferTask
  auto task = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  task.settings.states = {2, 2};
  task.settings.couple = {{1, 2}};
  task.run();

  // get values calculated in task
  auto calculator = task.getFDEETCalculator();
  auto H = calculator->getHamiltonian();
  auto S = calculator->getDeterminants();
  auto C = calculator->getLinCoeffs();
  auto E = calculator->getEigenValues();

  // check FDE-ET values
  EXPECT_NEAR(0.232994, calculator->getAnalyticalCoupling(), 1e-5); // delta E
  EXPECT_NEAR(0.528114, calculator->getLongRangeExEnergy(), 1e-5);  // V
  EXPECT_NEAR(-15.48758, H(0, 0), 1e-4);                            // H11
  EXPECT_NEAR(-15.31685, H(1, 1), 1e-4);                            // H22
  EXPECT_NEAR(-15.55357, H(0, 1), 1e-4);                            // H12
  EXPECT_NEAR(1.05024, S(0, 0), 1e-4);                              // S11
  EXPECT_NEAR(1.04470, S(1, 1), 1e-4);                              // S22
  EXPECT_NEAR(0.76112, S(0, 1), 1e-4);                              // S12
  EXPECT_NEAR(-0.82982, C(0, 0), 1e-4);                             // C11
  EXPECT_NEAR(1.15548, C(0, 1), 1e-4);                              // C12
  EXPECT_NEAR(-0.18945, C(1, 0), 1e-4);                             // C21
  EXPECT_NEAR(-1.41074, C(1, 1), 1e-4);                             // C22
  EXPECT_NEAR(-15.49697, E(0), 1e-4);                               // E1
  EXPECT_NEAR(-14.96886, E(1), 1e-4);                               // E2

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys1a->getSystemPath() + sys1a->getSystemName() + "_tmp/",
                                                        sys1a->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys2a->getSystemPath() + sys2a->getSystemName() + "_tmp/",
                                                        sys2a->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys1b->getSystemPath() + sys1b->getSystemName() + "_tmp/",
                                                        sys1b->getSystemName() + "_tmp");
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys2b->getSystemPath() + sys2b->getSystemName() + "_tmp/",
                                                        sys2b->getSystemName() + "_tmp");
  std::string basePath = sys1a->getSystemPath().substr(0, sys1a->getSystemPath().size() - sys1a->getSystemName().size() - 1);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(basePath + "sys-superSystem12/", "superSystem12");
}
/**
 * @test
 * @brief Tests disjoint functionality of ElectronTransferTask.h/.cpp.
 *        For the [BeH2]+ should yield exactly the same result, thus the same results are expected.
 *        Additionally run the printing of density matrix contributions just to check if the program does not die.
 */
TEST_F(ElectronTransferTaskTest, unrestrictedDisjoint) {
  // get Systems
  auto sys1a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  auto sys2a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  auto sys1b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  auto sys2b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // perform FeezeAndThawTask
  // state1
  auto task1 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1a, sys2a}, {});
  task1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task1.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task1.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task1.settings.maxCycles = 3;
  task1.run();

  // state2
  auto task2 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1b, sys2b}, {});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task2.settings.maxCycles = 3;
  task2.run();

  // perform ElectronTransferTask
  auto task = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  task.settings.states = {2, 2};
  task.settings.couple = {{1, 2}};
  task.settings.disjoint = {1, 2};
  task.settings.printContributions = true;
  task.run();

  // get values calculated in task
  auto calculator = task.getFDEETCalculator();
  auto H = calculator->getHamiltonian();
  auto S = calculator->getDeterminants();
  auto C = calculator->getLinCoeffs();
  auto E = calculator->getEigenValues();

  // check FDE-ET values
  EXPECT_NEAR(0.2205875, calculator->getAnalyticalCoupling(), 1e-5); // delta E
  EXPECT_NEAR(0.5058070, calculator->getLongRangeExEnergy(), 1e-5);  // V
  EXPECT_NEAR(-15.46291, H(0, 0), 1e-4);                             // H11
  EXPECT_NEAR(-15.29135, H(1, 1), 1e-4);                             // H22
  EXPECT_NEAR(-15.52435, H(0, 1), 1e-4);                             // H12
  EXPECT_NEAR(0.90410, S(0, 0), 1e-4);                               // S11
  EXPECT_NEAR(0.88999, S(1, 1), 1e-4);                               // S22
  EXPECT_NEAR(0.64632, S(0, 1), 1e-4);                               // S12
  EXPECT_NEAR(-0.90515, C(0, 0), 1e-4);                              // C11
  EXPECT_NEAR(1.22585, C(0, 1), 1e-4);                               // C12
  EXPECT_NEAR(-0.19320, C(1, 0), 1e-4);                              // C21
  EXPECT_NEAR(-1.51437, C(1, 1), 1e-4);                              // C22
  EXPECT_NEAR(-15.47110, E(0), 1e-4);                                // E1
  EXPECT_NEAR(-14.96529, E(1), 1e-4);                                // E2

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // clean up cubes
  std::vector<std::string> spin = {"alpha", "beta"};
  std::string basePath = sys1a->getSystemPath().substr(0, sys1a->getSystemPath().size() - sys1a->getSystemName().size() - 1);
  for (unsigned iState = 0; iState < 2; ++iState) {
    for (unsigned jState = 0; jState < 2; ++jState) {
      for (unsigned iSpin = 0; iSpin < 2; ++iSpin) {
        std::remove((basePath + "sys-superSystem-disjoint-12/diabDens" + std::to_string(iState) +
                     std::to_string(jState) + +"_" + spin[iSpin] + ".cube")
                        .c_str());
      }
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(basePath + "sys-superSystem-disjoint-12/",
                                                        "superSystem-disjoint-12");
}

/**
 * @test
 * @brief Tests the diskmode funtionality of ElectronTransferTask.h/.cpp.
 */
TEST_F(ElectronTransferTaskTest, diskMode) {
  // get Systems
  auto sys1a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  auto sys2a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  auto sys1b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  auto sys2b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // perform FeezeAndThawTask
  // state1
  auto task1 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1a, sys2a}, {});
  task1.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task1.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task1.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task1.settings.maxCycles = 3;
  task1.run();

  // state2
  auto task2 = FreezeAndThawTask<Options::SCF_MODES::UNRESTRICTED>({sys1b, sys2b}, {});
  task2.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  task2.settings.embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
  task2.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
  task2.settings.maxCycles = 3;
  task2.run();

  // perform ElectronTransferTask without diskmode
  auto task = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  task.settings.states = {2, 2};
  task.settings.couple = {{1, 2}};
  task.run();

  // perform ElectronTransferTask with diskmode
  auto taskDisk = ElectronTransferTask({sys1a, sys2a, sys1b, sys2b});
  taskDisk.settings.states = {2, 2};
  taskDisk.settings.couple = {{1, 2}};
  taskDisk.settings.diskMode = true;
  taskDisk.run();

  // get values calculated in tasks
  auto calculator = task.getFDEETCalculator();
  auto H = calculator->getHamiltonian();
  auto S = calculator->getDeterminants();
  auto C = calculator->getLinCoeffs();
  auto E = calculator->getEigenValues();

  auto calculatorDisk = taskDisk.getFDEETCalculator();
  auto HDisk = calculatorDisk->getHamiltonian();
  auto SDisk = calculatorDisk->getDeterminants();
  auto CDisk = calculatorDisk->getLinCoeffs();
  auto EDisk = calculatorDisk->getEigenValues();

  // check FDE-ET values
  EXPECT_NEAR(calculatorDisk->getAnalyticalCoupling(), calculator->getAnalyticalCoupling(), 1e-8); // delta E
  EXPECT_NEAR(calculatorDisk->getLongRangeExEnergy(), calculator->getLongRangeExEnergy(), 1e-8);   // V
  EXPECT_NEAR(HDisk(0, 0), H(0, 0), 1e-8);                                                         // H11
  EXPECT_NEAR(HDisk(1, 1), H(1, 1), 1e-8);                                                         // H22
  EXPECT_NEAR(HDisk(0, 1), H(0, 1), 1e-8);                                                         // H12
  EXPECT_NEAR(SDisk(0, 0), S(0, 0), 1e-8);                                                         // S11
  EXPECT_NEAR(SDisk(1, 1), S(1, 1), 1e-8);                                                         // S22
  EXPECT_NEAR(SDisk(0, 1), S(0, 1), 1e-8);                                                         // S12
  EXPECT_NEAR(CDisk(0, 0), C(0, 0), 1e-8);                                                         // C11
  EXPECT_NEAR(CDisk(0, 1), C(0, 1), 1e-8);                                                         // C12
  EXPECT_NEAR(CDisk(1, 0), C(1, 0), 1e-8);                                                         // C21
  EXPECT_NEAR(CDisk(1, 1), C(1, 1), 1e-8);                                                         // C22
  EXPECT_NEAR(EDisk(0), E(0), 1e-8);                                                               // E1
  EXPECT_NEAR(EDisk(1), E(1), 1e-8);                                                               // E2

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_1);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::BE_DEF2_TZVP_PW91_UNRES_0);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PW91_UNRES_1);

  // clean up t_dmats
  std::string basePath = sys1a->getSystemPath().substr(0, sys1a->getSystemPath().size() - sys1a->getSystemName().size() - 1);
  for (unsigned iState = 0; iState < 2; ++iState) {
    for (unsigned jState = 0; jState < 2; ++jState) {
      std::remove(
          (basePath + "sys-superSystem12/t_dmat" + std::to_string(iState) + std::to_string(jState) + ".dmat.unres.h5").c_str());
    }
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(basePath + "sys-superSystem12/", "superSystem12");
}
} /*namespace Serenity*/
