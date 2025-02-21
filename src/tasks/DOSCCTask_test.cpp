/**
 * @file DOSCCTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jul 9, 2021
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
#include "tasks/DOSCCTask.h"                          //To be tested.
#include "data/ElectronicStructure.h"                 //Get energy for comparison.
#include "io/FormattedOutputStream.h"                 //Change print level.
#include "math/FloatMaths.h"                          //Precision.
#include "settings/Settings.h"                        //Change system settings.
#include "system/SystemController.h"                  //Access to system settings.
#include "tasks/LocalCorrelationTask.h"               //DLPNO-CCSD(T0) reference calculations.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Access to test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class DOSCCTaskTest : public ::testing::Test {
 protected:
  DOSCCTaskTest() = default;

  virtual ~DOSCCTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
  void cleanUpA(std::shared_ptr<SystemController> system) {
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "a0/", "a0");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "a1/", "a1");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "a2/", "a2");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "a3/", "a3");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  }
  void cleanUpB(std::shared_ptr<SystemController> system) {
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "b0/", "b0");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "b1/", "b1");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "b2/", "b2");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "b3/", "b3");
    SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  }
};

TEST_F(DOSCCTaskTest, ethane_bond_stretching) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings aSettings = a->getSettings();
  aSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  aSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, aSettings);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  Settings bSettings = b->getSettings();
  bSettings.name = "b";
  bSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, bSettings);

  DOSCCTask dosccTask({a, b});
  dosccTask.settings.dosSettings = Options::DOS_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.printGroupAnalysis = true;
  dosccTask.settings.orbitalPairAnalysis = true;
  dosccTask.settings.lcPairSelection[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.strictTriples = false;
  dosccTask.run();
  const auto relativeEnergies = dosccTask.getRelativeEnergies();
  const auto totalEnergies = dosccTask.getTotalEnergies();

  // Small changes in the triple correction because of rigorous removal of distant orbital pairs.
  EXPECT_NEAR(relativeEnergies(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergies(1), 0.063639186813460924, NORMAL_D);
  EXPECT_NEAR(totalEnergies(0), -79.51949931969429, NORMAL_D);
  EXPECT_NEAR(totalEnergies(1), -79.45586013320775, NORMAL_D);

  const auto relativeEnergiesFromPairSelected = dosccTask.getRelativeEnergiesFromPairSelected();
  const auto totalEnergiesFromPairSelected = dosccTask.getTotalEnergiesFromPairSelected();

  EXPECT_NEAR(relativeEnergiesFromPairSelected(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergiesFromPairSelected(1), 0.063547338789121, NORMAL_D);
  EXPECT_NEAR(totalEnergiesFromPairSelected(0), -79.521960798079562, NORMAL_D);
  EXPECT_NEAR(totalEnergiesFromPairSelected(1), -79.45841346071019, NORMAL_D);

  cleanUpA(a);
  cleanUpB(b);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DOSCCTaskTest, ethane_bond_stretching_vsFullCC_pairSelected) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings aSettings = a->getSettings();
  aSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  aSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, aSettings);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  Settings bSettings = b->getSettings();
  bSettings.name = "b";
  bSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, bSettings);

  DOSCCTask dosccTask({a, b});
  dosccTask.settings.dosSettings = Options::DOS_SETTINGS::LOOSE;
  dosccTask.settings.alignOrbitals = false;
  dosccTask.settings.wfemb.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[0].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[1].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[2].pnoSettings = Options::PNO_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[2].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[3].pnoSettings = Options::PNO_SETTINGS::LOOSE;
  dosccTask.settings.wfemb.lcSettings[3].useFrozenCore = true;
  dosccTask.settings.printGroupAnalysis = true;
  dosccTask.settings.orbitalPairAnalysis = true;
  dosccTask.settings.lcPairSelection[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[0].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[0].enforceHFFockian = true;
  dosccTask.settings.lcPairSelection[1].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[1].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[2].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[2].useFrozenCore = true;
  dosccTask.run();

  const auto relativeEnergiesFromPairSelected = dosccTask.getRelativeEnergiesFromPairSelected();
  const auto totalEnergiesFromPairSelected = dosccTask.getTotalEnergiesFromPairSelected();
  LocalCorrelationTask lcTaskA(a);
  lcTaskA.settings.lcSettings.useFrozenCore = true;
  lcTaskA.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskA.run();
  const double energyA = a->getElectronicStructure<RESTRICTED>()->getEnergy();
  LocalCorrelationTask lcTaskB(b);
  lcTaskB.settings.lcSettings.useFrozenCore = true;
  lcTaskB.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskB.run();
  const double energyB = b->getElectronicStructure<RESTRICTED>()->getEnergy();
  const double relativeEnergy = energyB - energyA;

  EXPECT_NEAR(relativeEnergiesFromPairSelected(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergiesFromPairSelected(1), relativeEnergy, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(0), energyA, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(1), energyB, 1e-8);

  cleanUpA(a);
  cleanUpB(b);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DOSCCTaskTest, ethane_bond_stretching_vsFullCC_pairSelected_lmp2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings aSettings = a->getSettings();
  aSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  aSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, aSettings);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  Settings bSettings = b->getSettings();
  bSettings.name = "b";
  bSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, bSettings);

  DOSCCTask dosccTask({a, b});
  dosccTask.settings.dosSettings = Options::DOS_SETTINGS::SPREAD;
  dosccTask.settings.alignOrbitals = false;
  dosccTask.settings.wfemb.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[0].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[0].method = Options::PNO_METHOD::DLPNO_MP2;
  dosccTask.settings.wfemb.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[1].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[1].method = Options::PNO_METHOD::DLPNO_MP2;
  dosccTask.settings.wfemb.lcSettings[2].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[2].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[2].method = Options::PNO_METHOD::DLPNO_MP2;
  dosccTask.settings.wfemb.lcSettings[3].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[3].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[3].method = Options::PNO_METHOD::DLPNO_MP2;
  dosccTask.settings.printGroupAnalysis = true;
  dosccTask.settings.orbitalPairAnalysis = true;
  dosccTask.settings.lcPairSelection[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[0].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[0].enforceHFFockian = true;
  dosccTask.settings.lcPairSelection[1].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[1].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[2].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[2].useFrozenCore = true;
  dosccTask.run();

  const auto relativeEnergiesFromPairSelected = dosccTask.getRelativeEnergiesFromPairSelected();
  const auto totalEnergiesFromPairSelected = dosccTask.getTotalEnergiesFromPairSelected();
  LocalCorrelationTask lcTaskA(a);
  lcTaskA.settings.lcSettings.useFrozenCore = true;
  lcTaskA.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskA.run();
  const double energyA = a->getElectronicStructure<RESTRICTED>()->getEnergy();
  LocalCorrelationTask lcTaskB(b);
  lcTaskB.settings.lcSettings.useFrozenCore = true;
  lcTaskB.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskB.run();
  const double energyB = b->getElectronicStructure<RESTRICTED>()->getEnergy();
  const double relativeEnergy = energyB - energyA;

  EXPECT_NEAR(relativeEnergiesFromPairSelected(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergiesFromPairSelected(1), relativeEnergy, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(0), energyA, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(1), energyB, 1e-8);

  cleanUpA(a);
  cleanUpB(b);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DOSCCTaskTest, ethane_bond_stretching_vsFullCC_pairSelected_scmp2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings aSettings = a->getSettings();
  aSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  aSettings.name = "a";
  a = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, aSettings);
  auto b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, true);
  Settings bSettings = b->getSettings();
  bSettings.name = "b";
  bSettings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  b = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneB_Def2_SVP_BP86, bSettings);

  DOSCCTask dosccTask({a, b});
  dosccTask.settings.dosSettings = Options::DOS_SETTINGS::SPREAD;
  dosccTask.settings.alignOrbitals = false;
  dosccTask.settings.wfemb.lcSettings[0].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[0].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[0].method = Options::PNO_METHOD::SC_MP2;
  dosccTask.settings.wfemb.lcSettings[1].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[1].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[1].method = Options::PNO_METHOD::SC_MP2;
  dosccTask.settings.wfemb.lcSettings[2].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[2].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[2].method = Options::PNO_METHOD::SC_MP2;
  dosccTask.settings.wfemb.lcSettings[3].pnoSettings = Options::PNO_SETTINGS::TIGHT;
  dosccTask.settings.wfemb.lcSettings[3].useFrozenCore = true;
  dosccTask.settings.wfemb.lcSettings[3].method = Options::PNO_METHOD::SC_MP2;
  dosccTask.settings.printGroupAnalysis = true;
  dosccTask.settings.orbitalPairAnalysis = true;
  dosccTask.settings.lcPairSelection[0].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[0].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[0].enforceHFFockian = true;
  dosccTask.settings.lcPairSelection[1].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[1].useFrozenCore = true;
  dosccTask.settings.lcPairSelection[2].method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  dosccTask.settings.lcPairSelection[2].useFrozenCore = true;
  dosccTask.settings.skipCrudePresPairSelected = true;
  dosccTask.run();

  const auto relativeEnergiesFromPairSelected = dosccTask.getRelativeEnergiesFromPairSelected();
  const auto totalEnergiesFromPairSelected = dosccTask.getTotalEnergiesFromPairSelected();
  LocalCorrelationTask lcTaskA(a);
  lcTaskA.settings.lcSettings.useFrozenCore = true;
  lcTaskA.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskA.run();
  const double energyA = a->getElectronicStructure<RESTRICTED>()->getEnergy();
  LocalCorrelationTask lcTaskB(b);
  lcTaskB.settings.lcSettings.useFrozenCore = true;
  lcTaskB.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD_T0;
  lcTaskB.run();
  const double energyB = b->getElectronicStructure<RESTRICTED>()->getEnergy();
  const double relativeEnergy = energyB - energyA;

  EXPECT_NEAR(relativeEnergiesFromPairSelected(0), 0.0, TIGHT_D);
  EXPECT_NEAR(relativeEnergiesFromPairSelected(1), relativeEnergy, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(0), energyA, 1e-8);
  EXPECT_NEAR(totalEnergiesFromPairSelected(1), energyB, 1e-8);

  cleanUpA(a);
  cleanUpB(b);
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

} /* namespace Serenity */
