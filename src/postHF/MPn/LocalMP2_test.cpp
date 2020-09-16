/**
 * @file LocalMP2_test.cpp
 *
 * @date Dec 14, 2018
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

#ifndef POSTHF_MPN_LOCALMP2_TEST_CPP_
#define POSTHF_MPN_LOCALMP2_TEST_CPP_

/* Include Serenity Internal Headers */
#include "postHF/MPn/LocalMP2.h"                      //To be tested.
#include "data/ElectronicStructure.h"                 //Run scf.
#include "postHF/MPn/MP2.h"                           //Reference.
#include "postHF/MPn/RIMP2.h"                         //Reference.
#include "system/SystemController.h"                  //Test systems.
#include "tasks/LocalizationTask.h"                   //Orbital localization.
#include "tasks/TDEmbeddingTask.h"                    //Embedded LMP2.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class LocalMP2Test : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Minimal test example. Only one orbital pair.
 */
TEST_F(LocalMP2Test, H2) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;
  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.02;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-5;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  LocalMP2 localMP2Calculator(lCController);
  RIMP2<scfMode> canonicalMP2(system);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  double canonicalMP2Energy = canonicalMP2.calculateCorrection();
  EXPECT_NEAR(canonicalMP2Energy, localMP2Energy, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test
 * @brief Minimal test example. Only one orbital pair.
 */
TEST_F(LocalMP2Test, H2FullFourCenter) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;
  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.02;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-5;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  LocalMP2 localMP2Calculator(lCController);
  localMP2Calculator.settings.useFourCenterIntegrals = true;
  MP2EnergyCorrector<scfMode> canonicalMP2(system);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  double canonicalMP2Energy = canonicalMP2.calculateElectronicEnergy();
  EXPECT_NEAR(canonicalMP2Energy, localMP2Energy, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test
 * @brief Multiple orbital pairs.
 */
TEST_F(LocalMP2Test, Water) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;
  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-5;
  lCSettings.pnoThreshold = 1e-12;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  LocalMP2 localMP2Calculator(lCController);

  RIMP2<scfMode> canonicalMP2(system);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  double canonicalMP2Energy = canonicalMP2.calculateCorrection();
  EXPECT_NEAR(canonicalMP2Energy, localMP2Energy, 1e-9);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
}

/**
 * @test
 * @brief Multiple orbital pairs.
 */
TEST_F(LocalMP2Test, WaterReduced_PAO_space) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;

  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.02;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-3;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  LocalMP2 localMP2Calculator(lCController);

  RIMP2<scfMode> canonicalMP2(system);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  double canonicalMP2Energy = canonicalMP2.calculateCorrection();
  EXPECT_NEAR(canonicalMP2Energy, localMP2Energy, 5e-4);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
}

/**
 * @test
 * @brief Multiple orbital pairs and localized orbitals.
 */
TEST_F(LocalMP2Test, LocalizedWater) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;
  RIMP2<scfMode> canonicalMP2(system);
  double canonicalMP2Energy = canonicalMP2.calculateCorrection();
  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-5;
  lCSettings.fockMatrixPreescreeningThreshold = 1e-5;
  lCSettings.pnoThreshold = 1e-12;
  lCSettings.doiPAOThreshold = 1e-7;
  lCSettings.diisStartResidual = 1e-4;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  LocalMP2 localMP2Calculator(lCController);
  localMP2Calculator.settings.maxResidual = 1e-5;

  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  EXPECT_NEAR(canonicalMP2Energy, localMP2Energy, 1e-4);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
}

/**
 * @test
 * @brief Embedding example in supersystem basis.
 */
TEST_F(LocalMP2Test, WaterDimer) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);
  TDEmbeddingTask<scfMode> tdTask(act, env);
  tdTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  tdTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  tdTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  lCSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  std::vector<std::shared_ptr<SystemController>> tmp = {env};
  auto lCController = std::make_shared<LocalCorrelationController>(act, lCSettings, tmp);
  LocalMP2 localMP2Calculator(lCController);

  RIMP2<scfMode> canonicalMP2(act);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  EXPECT_NEAR(-0.20587434893492509, localMP2Energy, 1e-5);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
}

/**
 * @test
 * @brief Embedding example in truncated basis.
 */
TEST_F(LocalMP2Test, WaterDimer_BasisTruncation) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP, true);
  TDEmbeddingTask<scfMode> tdTask(act, env);
  tdTask.settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  tdTask.settings.embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  tdTask.settings.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION;
  tdTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.embeddingSettings.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PBE;
  lCSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::LEVELSHIFT;
  std::vector<std::shared_ptr<SystemController>> tmp = {env};
  auto lCController = std::make_shared<LocalCorrelationController>(act, lCSettings, tmp);
  LocalMP2 localMP2Calculator(lCController);

  RIMP2<scfMode> canonicalMP2(act);
  double localMP2Energy = localMP2Calculator.calculateEnergyCorrection().sum();
  EXPECT_NEAR(-0.2060700782541677, localMP2Energy, 5e-4);

  std::string name = act->getSystemName() + "+" + env->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(env->getSystemPath() + name + "/", name);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
}

} /* namespace Serenity */

#endif /* POSTHF_MPN_LOCALMP2_TEST_CPP_ */
