/**
 * @file DLPNO_CCSD_test.cpp
 *
 * @date Jun 25, 2019
 * @author Moritz Bensberg
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
#include "postHF/CC/DLPNO_CCSD.h" //To be tested.
#include "data/ElectronicStructure.h"
#include "energies/EnergyContributions.h"
#include "postHF/CC/CCSD_T.h"                                   //Reference.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Access to local-correlation framework.
#include "settings/Settings.h"                                  //Test systems.
#include "system/SystemController.h"                            //Test systems.
#include "tasks/CoupledClusterTask.h"                           //Test via CC-task.
#include "tasks/DFTEmbeddedLocalCorrelationTask.h"              //Test for embedded CCSD calculation.
#include "tasks/LocalizationTask.h"                             //Orbital localization.
#include "tasks/ScfTask.h"                                      //Clean electronic structures.
#include "testsupply/SystemController__TEST_SUPPLY.h"           //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <iomanip>       //std::setprecision

namespace Serenity {
class DLPNO_CCSDTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(DLPNO_CCSDTest, H2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  auto settings = system->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.dumpIntegrals = false;
  lCSettings.diisStartResidual = 1e+100;
  lCSettings.dampingFactor = 0.0;
  lCSettings.finalDamping = 0.0;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-9, 100);
  localCCSD.testRun = true;
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 2.4e-9, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, H2_RI) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-7, 100);
  //  localCCSD.testRun = true;
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-9, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, H2Dimer) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP, true);
  auto settings = system->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP, settings);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.dumpIntegrals = false;
  lCSettings.diisStartResidual = 1e+100;
  lCSettings.dampingFactor = 0.0;
  lCSettings.finalDamping = 0.0;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 5e-10, 100);
  localCCSD.testRun = true;
  std::cout << std::setprecision(10);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-7, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, H2Dimer_RI) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 0.0;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-5, 100);
  //  localCCSD.testRun = true;
  //  localCCSD.invertTest = true;
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-5, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, WaterFull4Center) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  Settings settings = system->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 0.0;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 0.0;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-6, 100);
  localCCSD.testRun = true;
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-5, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, Water_RI) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  Settings settings = system->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, settings);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 0.0;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 4e-5, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-5, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy.second << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy.second << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy.second, 3e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, Water_RI_Localized_IBO) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-5, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections().second;
  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 3.33e-7;
  lCSettings.doiPAOThreshold = 1e-2;
  lCSettings.mullikenThreshold = 1e-3;
  lCSettings.orbitalToShellThreshold = 1e-3;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-7, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 2e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, Water_RI_Localized_FB) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  // This object can do both, CCSD(T) and CCSD.
  CCSD_T canonicalCC(system, 1e-5, 100);
  auto canonicalCCSDEnergy = canonicalCC.calculateElectronicEnergyCorrections().second;
  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 3.33e-7;
  lCSettings.doiPAOThreshold = 1e-2;
  lCSettings.mullikenThreshold = 1e-3;
  lCSettings.orbitalToShellThreshold = 1e-3;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-7, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 2e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, Ethane) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.dumpIntegrals = false;
  lCSettings.diisStartResidual = 1e-10;
  lCSettings.dampingFactor = 0.0;
  lCSettings.finalDamping = 0.0;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-7, 100);
  localCCSD.testRun = true;
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  // Canonical reference is taken from Serenity. Calculating the reference on the fly is to slow
  // for my liking.
  double canonicalCCSDEnergy = -3.456010942618e-01;
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, Ethane_Localized) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  // Canonical reference is taken from Serenity. Calculating the reference on the fly is to slow
  // for my liking.
  double canonicalCCSDEnergy = -3.456010942618e-01;
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-7;
  lCSettings.mullikenThreshold = 1e-3;
  lCSettings.orbitalToShellThreshold = 1e-3;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, Ethane_Localized_linSig) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  // Canonical reference is taken from Serenity. Calculating the reference on the fly is to slow
  // for my liking.
  double canonicalCCSDEnergy = -3.456010942618e-01;
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
  lCSettings.dumpIntegrals = false;
  lCSettings.linearScalingSigmaVector = true;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, Ethane_Localized_dumpIntegrals) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, true);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  // Canonical reference is taken from Serenity. Calculating the reference on the fly is to slow
  // for my liking.
  double canonicalCCSDEnergy = -3.456010942618e-01;
  LocalizationTask locTask(system);
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-7;
  lCSettings.mullikenThreshold = 1e-3;
  lCSettings.orbitalToShellThreshold = 1e-3;
  lCSettings.dumpIntegrals = true;
  lCSettings.maximumMemoryRatio = 1e-3;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  double ccsdEnergy = localCCSD.calculateElectronicEnergyCorrections().sum();
  std::cout << "Canonical CCSD Energy: " << canonicalCCSDEnergy << std::endl;
  std::cout << "Delta E: " << ccsdEnergy - canonicalCCSDEnergy << std::endl;
  EXPECT_NEAR(0.0, ccsdEnergy - canonicalCCSDEnergy, 5e-4);
  auto orbitalPairs = localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, PhenolEmbedded) {
  SystemController__TEST_SUPPLY::cleanUp();
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  std::vector<std::shared_ptr<SystemController>> tmp = {env};
  DFTEmbeddedLocalCorrelationTask dftEmbeddedLocalCorrelationTask(act, tmp);
  dftEmbeddedLocalCorrelationTask.settings.fromSupersystem = true;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  dftEmbeddedLocalCorrelationTask.settings.trunc.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  dftEmbeddedLocalCorrelationTask.settings.trunc.netThreshold = 1e-4;
  dftEmbeddedLocalCorrelationTask.settings.normThreshold = 1e-5;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD;
  dftEmbeddedLocalCorrelationTask.settings.split.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD;
  dftEmbeddedLocalCorrelationTask.run();
  double ccsdEnergy = act->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  EXPECT_NEAR(ccsdEnergy, -0.21755530254232638, 1e-5);
  std::string name = "TMP_Supersystem";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, PhenolEmbedded_dftOrbitals) {
  SystemController__TEST_SUPPLY::cleanUp();
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  act->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::DFT);
  env->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::DFT);
  std::vector<std::shared_ptr<SystemController>> tmp = {env};
  DFTEmbeddedLocalCorrelationTask dftEmbeddedLocalCorrelationTask(act, tmp);
  dftEmbeddedLocalCorrelationTask.settings.fromSupersystem = true;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.enforceHFFockian = true;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.useProjectedOccupiedOrbitals = true;
  dftEmbeddedLocalCorrelationTask.settings.normThreshold = 1e-5;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD;
  dftEmbeddedLocalCorrelationTask.settings.split.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD;
  dftEmbeddedLocalCorrelationTask.run();
  double ccsdEnergy = act->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  EXPECT_NEAR(ccsdEnergy, -0.22491401778086129, 1e-5);
  std::string name = "TMP_Supersystem";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, PhenolEmbedded_projectedOccupiedOrbitals) {
  SystemController__TEST_SUPPLY::cleanUp();
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  act->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::DFT);
  env->setElectronicStructureMethod(Options::ELECTRONIC_STRUCTURE_THEORIES::DFT);
  std::vector<std::shared_ptr<SystemController>> tmp = {env};
  DFTEmbeddedLocalCorrelationTask dftEmbeddedLocalCorrelationTask(act, tmp);
  dftEmbeddedLocalCorrelationTask.settings.fromSupersystem = true;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.enforceHFFockian = true;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::HUZINAGA;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  dftEmbeddedLocalCorrelationTask.settings.normThreshold = 1e-5;
  dftEmbeddedLocalCorrelationTask.settings.lcSettings.method = Options::PNO_METHOD::DLPNO_CCSD;
  dftEmbeddedLocalCorrelationTask.settings.split.systemPartitioning = Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD;
  dftEmbeddedLocalCorrelationTask.run();
  double ccsdEnergy = act->getElectronicStructure<RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::CCSD_CORRECTION);
  EXPECT_NEAR(ccsdEnergy, -0.2202521621686124, 1e-5);
  std::string name = "TMP_Supersystem";
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(DLPNO_CCSDTest, Phenol) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  auto phenol = *act + *env;
  ScfTask<RESTRICTED> scfTask(phenol);
  scfTask.run();
  // LOOSE PNO
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  lCSettings.linearScalingSigmaVector = false;

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;

  LocalizationTask locTask(phenol);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  CoupledClusterTask ccTask(phenol);
  ccTask.settings.level = Options::CC_LEVEL::DLPNO_CCSD;
  ccTask.settings.lcSettings = lCSettings;
  ccTask.run();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  // Canonical result ca. -306.3553068128
  EXPECT_NEAR(phenol->getElectronicStructure<RESTRICTED>()->getEnergy(), -306.36245921105103, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(phenol);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DLPNO_CCSDTest, Phenol_linSig) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  auto phenol = *act + *env;
  ScfTask<RESTRICTED> scfTask(phenol);
  scfTask.run();
  // LOOSE PNO
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  lCSettings.linearScalingSigmaVector = true;

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;

  LocalizationTask locTask(phenol);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();
  CoupledClusterTask ccTask(phenol);
  ccTask.settings.level = Options::CC_LEVEL::DLPNO_CCSD;
  ccTask.settings.lcSettings = lCSettings;
  ccTask.run();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  // Canonical result ca. -306.3553068128
  EXPECT_NEAR(phenol->getElectronicStructure<RESTRICTED>()->getEnergy(), -306.362790, 5e-6);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(phenol);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
