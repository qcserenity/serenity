/**
 * @file DLPNO_CCSD_T0_test.cpp
 *
 * @author Moritz Bensberg
 * @date Nov 8, 2019
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
#include "postHF/CC/DLPNO_CCSD_T0.h"                            //To be tested.
#include "io/FormattedOutputStream.h"                           //Print levels.
#include "postHF/CC/CCSD_T.h"                                   //References.
#include "postHF/CC/DLPNO_CCSD.h"                               //Initial CCSD calculation.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Easy access to local-correlation framework.
#include "system/SystemController.h"                            //Test systems.
#include "tasks/LocalizationTask.h"                             //Orbital localization.
#include "tasks/ScfTask.h"                                      //Clean electronic structures.
#include "testsupply/SystemController__TEST_SUPPLY.h"           //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class DLPNO_CCSD_T0Test : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests the DLPNO-CCSD(T0) for water and canonic orbitals. Canonical CCSD(T) is used as a reference.
 */
TEST_F(DLPNO_CCSD_T0Test, WaterCanonic) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 8e-4, 100);
  localCCSD.dlpnoCCSDSettings.keepMO3CenterIntegrals = true;
  double ccsdCorrection = localCCSD.calculateElectronicEnergyCorrections().sum();
  double triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);

  // Reference
  CCSD_T canonicalCCSD_T(system, 1e-5, 100);
  double canonicalCCSD = canonicalCCSD_T.calculateElectronicEnergyCorrections().second;
  double canonicalTriples = canonicalCCSD_T.calculateTripplesCorrection();

  EXPECT_NEAR(0.0, canonicalCCSD - ccsdCorrection, 3e-4);
  EXPECT_NEAR(0.0, canonicalTriples - triplesCorrection, 1e-4);

  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the DLPNO-CCSD(T0) for water and local orbitals.
 */
TEST_F(DLPNO_CCSD_T0Test, WaterLocal) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();

  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  // Local CCSD_T0
  LocalCorrelationSettings lCSettings;
  lCSettings.diisStartResidual = 1.0; // Convergence is quite hard for water ...
  lCSettings.pnoThreshold = 1e-7;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.tnoThreshold = 1e-9;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  localCCSD.dlpnoCCSDSettings.keepMO3CenterIntegrals = true;
  double ccsdCorrection = localCCSD.calculateElectronicEnergyCorrections().sum();
  double triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
  OutputControl::vOut << "CCSD " << ccsdCorrection << std::endl;
  OutputControl::vOut << "T0   " << triplesCorrection << std::endl;
  EXPECT_NEAR(-0.213858, ccsdCorrection, 1e-4);
  EXPECT_NEAR(-0.0029991450991246017, triplesCorrection, 1e-4);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}
/**
 * @test
 * @brief Tests the DLPNO-CCSD(T0) for water and local orbitals with loose thresholds.
 *        Check the triples lists explicitly.
 */
TEST_F(DLPNO_CCSD_T0Test, WaterLocal_LOOSE) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();

  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  // Local CCSD_T0
  LocalCorrelationSettings lCSettings;
  lCSettings.tnoThreshold = 1e-9;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  localCCSD.dlpnoCCSDSettings.keepMO3CenterIntegrals = true;
  double ccsdCorrection = localCCSD.calculateElectronicEnergyCorrections().sum();
  double triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
  auto triples = localCorrelationController->getOrbitalTriples();
  auto distantTriplesPairs = localCorrelationController->getOrbitalPairs(OrbitalPairTypes::DISTANT_TRIPLES);
  auto closePairs = localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  OutputControl::vOut << "CCSD " << ccsdCorrection << std::endl;
  OutputControl::vOut << "T0   " << triplesCorrection << std::endl;
  EXPECT_NEAR(-0.21405024590820307, ccsdCorrection, 1e-4);
  EXPECT_NEAR(-0.0029991450991246017, triplesCorrection, 1e-4);
  EXPECT_EQ(triples.size(), 16);
  EXPECT_EQ(closePairs.size(), 11);
  EXPECT_EQ(distantTriplesPairs.size(), 4);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the DLPNO-CCSD(T0) for ethane and canonic orbitals. Canonical CCSD(T) is used as a reference.
 */
TEST_F(DLPNO_CCSD_T0Test, EthaneCanonic) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-5, 100);
  localCCSD.dlpnoCCSDSettings.keepMO3CenterIntegrals = true;
  double ccsdCorrection = localCCSD.calculateElectronicEnergyCorrections().sum();
  double triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);

  // Reference
  CCSD_T canonicalCCSD_T(system, 1e-5, 100);
  double canonicalCCSD = canonicalCCSD_T.calculateElectronicEnergyCorrections().second;
  double canonicalTriples = canonicalCCSD_T.calculateTripplesCorrection();

  EXPECT_NEAR(0.0, canonicalCCSD - ccsdCorrection, 3e-4);
  EXPECT_NEAR(0.0, canonicalTriples - triplesCorrection, 3e-4);

  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests the DLPNO-CCSD(T0) for ethane and local orbitals.
 *        In this example weak pairs occure, thus triples prescreening will be used.
 */
TEST_F(DLPNO_CCSD_T0Test, EthaneLocal) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();

  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  // Local CCSD_T0
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-5, 100);
  localCCSD.dlpnoCCSDSettings.keepMO3CenterIntegrals = true;
  double ccsdCorrection = localCCSD.calculateElectronicEnergyCorrections().sum();
  double triplesCorrection = DLPNO_CCSD_T0::calculateEnergyCorrection(localCorrelationController);
  OutputControl::vOut << "CCSD " << ccsdCorrection << std::endl;
  OutputControl::vOut << "T0   " << triplesCorrection << std::endl;
  EXPECT_NEAR(-0.34562194336802227, ccsdCorrection, 1e-5);
  EXPECT_NEAR(-0.0084892772234549819, triplesCorrection, 1e-5);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
