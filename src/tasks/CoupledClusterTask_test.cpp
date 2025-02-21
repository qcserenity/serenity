/**
 * @file CoupledClusterTask_test.cpp
 *
 * @date Nov 17, 2017
 * @author Jan Unsleber
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
#include "tasks/CoupledClusterTask.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <istream>

namespace Serenity {
/**
 * @class CoupledClusterTaskTest
 * @brief Sets everything up for the tests of CoupledClusterTaskTest.h/.cpp .
 */
class CoupledClusterTaskTest : public ::testing::Test {
 protected:
  CoupledClusterTaskTest() {
  }

  virtual ~CoupledClusterTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests restricted output.
 */
TEST_F(CoupledClusterTaskTest, minimal) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  CoupledClusterTask task(act);
  task.run();
}

TEST_F(CoupledClusterTaskTest, WaterDLPNOCCSDT0_TIGHT) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  iOOptions.timingsPrintLevel = 1;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalizationTask loc(act);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.settings.splitValenceAndCore = true;
  loc.run();
  CoupledClusterTask task(act);
  task.settings.level = Options::CC_LEVEL::DLPNO_CCSD_T0;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(-0.215243, ccsdCorrection, 1e-6);
  // Canonical: -0.00329571. However, this is only semi-canonical treatment!
  EXPECT_NEAR(-0.003178, triplesCorrection, 1e-6);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(CoupledClusterTaskTest, WaterDLPNOCCSDT0_NORMAL) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  iOOptions.timingsPrintLevel = 1;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalizationTask loc(act);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.settings.splitValenceAndCore = true;
  loc.run();
  CoupledClusterTask task(act);
  task.settings.level = Options::CC_LEVEL::DLPNO_CCSD_T0;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::NORMAL;
  task.settings.lcSettings.ignoreMemoryConstraints = true;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(-0.215187, ccsdCorrection, 1e-6);
  // Canonical: -0.00329571. However, this is only semi-canonical treatment!
  EXPECT_NEAR(-0.003182, triplesCorrection, 1e-6);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(CoupledClusterTaskTest, WaterDLPNOCCSDT0_LOOSE) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  iOOptions.timingsPrintLevel = 1;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<RESTRICTED> scf(act);
  scf.run();
  LocalizationTask loc(act);
  loc.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  loc.settings.splitValenceAndCore = true;
  loc.run();
  CoupledClusterTask task(act);
  task.settings.level = Options::CC_LEVEL::DLPNO_CCSD_T0;
  task.settings.normThreshold = 1e-6;
  task.settings.lcSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  task.run();
  double ccsdCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  double triplesCorrection =
      act->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(-0.215378, ccsdCorrection, 1e-6);
  // Canonical: -0.00329571. However, this is only semi-canonical treatment!
  EXPECT_NEAR(-0.0031526883, triplesCorrection, 1e-6);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
