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
  lCSettings.fockMatrixPrescreeningThresholdd = 1e-5;
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
  lCSettings.useProjectedOccupiedOrbitals = true;
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

/**
 * @test
 * @brief Unrelaxed MP2 density.
 */
TEST_F(LocalMP2Test, Water_Density) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  double e = system->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;

  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-5;
  lCSettings.pnoThreshold = 1e-12;

  auto localSystem = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  double eLocal = localSystem->getElectronicStructure<scfMode>()->getEnergy();
  (void)eLocal;
  auto lCController = std::make_shared<LocalCorrelationController>(localSystem, lCSettings);
  LocalMP2 localMP2Calculator(lCController);
  localMP2Calculator.calculateEnergyCorrection();
  auto localDensity = localMP2Calculator.calculateDensityCorrection();

  Eigen::MatrixXd refLocalDensity(24, 24);
  refLocalDensity << 0.0113596, 0.00076305, 1.49736e-11, -0.00123388, 0.00021755, 0.0016788, -0.000352098, 1.92926e-11,
      -0.000342357, -0.000642808, -0.00213377, -0.00608358, -0.00431697, -1.0781e-11, -0.0190573, 0.000138104,
      3.35572e-10, -0.00623952, -0.000190833, -4.00913e-11, 9.93093e-11, 0.00133348, -0.000374311, 0.000513156,
      0.00076305, 0.00218093, -9.10822e-11, 0.000663383, 0.00016383, -0.000352098, -0.00090616, -2.16365e-11,
      0.000174876, -0.000114754, 0.00021921, 0.00155907, -0.00231078, 1.2152e-09, -0.00675013, 0.000409342, -6.79241e-10,
      0.00114005, 0.000784914, 1.25614e-11, -3.2557e-11, -0.000219469, -2.23044e-05, 0.000184521, 1.49736e-11,
      -9.10822e-11, 0.00117035, -1.67114e-11, 5.66701e-12, 7.56952e-11, -1.03888e-10, -0.000113203, -5.35077e-12,
      -2.10583e-11, -5.28659e-12, -7.65893e-11, 1.53093e-10, -0.00167752, 5.50973e-11, -6.39245e-11, 0.000483885,
      -3.09858e-11, 6.48358e-11, -2.41836e-06, 0.000868814, -8.87137e-12, 4.90901e-12, -7.28254e-12, -0.00123388,
      0.000663383, -1.67114e-11, 0.000915778, 2.37111e-05, -0.000531802, -0.000156243, -7.38655e-12, 3.6619e-05,
      0.000145493, 0.000446853, 0.00170266, -0.000127843, 1.01179e-10, 0.00148885, -8.35569e-05, -1.21659e-10,
      0.00107889, 0.000498621, 6.05172e-12, -2.51755e-11, -0.000144674, 5.55284e-05, -7.92728e-05, 0.00021755,
      0.00016383, 5.66701e-12, 2.37111e-05, 0.00102696, -0.000497592, 0.000139061, 1.03535e-11, 0.000277902, 0.000150518,
      4.54019e-05, 0.000293763, -0.000313879, -9.08818e-11, -0.000705467, -0.00125455, 8.38774e-11, 5.94147e-06,
      0.000290727, 2.69224e-12, -3.07162e-13, 4.67149e-05, 0.000727652, 0.000115256, 0.0016788, -0.000352098,
      7.56952e-11, -0.000531802, -0.000497592, 0.0113596, 0.00076305, -1.27561e-11, 0.000530617, -0.00113501,
      -0.00213377, -0.00608357, -0.00431697, 1.32054e-09, 0.00508405, -0.0183672, 1.20741e-09, 0.00143661, -0.00607488,
      1.21818e-11, 2.38092e-11, 4.52972e-05, -0.000774462, 0.00125689, -0.000352098, -0.00090616, -1.03888e-10,
      -0.000156243, 0.000139061, 0.00076305, 0.00218093, -2.18383e-10, -1.41277e-05, 0.000683168, 0.000219209,
      0.00155907, -0.00231078, 4.95403e-09, 0.00214883, -0.00641206, -2.09349e-09, 0.000461806, 0.00130481, -6.3105e-11,
      3.83396e-11, 0.000246225, 0.000122354, -8.43477e-05, 1.92926e-11, -2.16365e-11, -0.000113203, -7.38655e-12,
      1.03535e-11, -1.27561e-11, -2.18383e-10, 0.00117034, 1.88905e-11, -5.26603e-11, 1.93433e-11, 5.30204e-12,
      1.48378e-10, -0.00167752, -3.95918e-11, 1.08879e-10, 0.000483885, -2.24067e-11, 9.61837e-11, 0.000838358,
      -0.000228034, 2.28275e-12, 4.40587e-13, 2.03255e-11, -0.000342357, 0.000174876, -5.35077e-12, 3.6619e-05,
      0.000277902, 0.000530617, -1.41277e-05, 1.88905e-11, 0.00100756, 4.84015e-05, -7.22395e-05, -0.000158636,
      -0.000269892, -9.16991e-11, -0.000871478, -0.00134043, 1.60723e-10, 0.000217342, -0.000225828, -3.90741e-12,
      -2.24048e-11, -0.000225506, 0.000592509, 0.000309832, -0.000642808, -0.000114754, -2.10583e-11, 0.000145493,
      0.000150518, -0.00113501, 0.000683168, -5.26603e-11, 4.84015e-05, 0.000935177, 0.000443306, 0.00172051,
      -0.000204992, 3.90207e-10, -0.000718517, 0.00110577, -4.27875e-10, 0.000266852, 0.00115228, -1.45357e-11,
      -2.77408e-12, -9.22014e-05, 0.000253637, -6.70315e-05, -0.00213377, 0.00021921, -5.28659e-12, 0.000446853,
      4.54019e-05, -0.00213377, 0.000219209, 1.93433e-11, -7.22395e-05, 0.000443306, 6.14049e-05, 0.00360498,
      -0.00530849, -3.87718e-10, 0.000948101, 0.00123686, -8.56294e-12, 0.000904666, 0.0011802, 2.3375e-11, -8.69035e-12,
      -5.15423e-05, 7.38745e-05, -4.95276e-05, -0.00608358, 0.00155907, -7.65893e-11, 0.00170266, 0.000293763,
      -0.00608357, 0.00155907, 5.30204e-12, -0.000158636, 0.00172051, 0.00360498, 0.0236029, -0.0252191, -8.02206e-10,
      0.00156873, 0.00204651, -8.7714e-10, 0.00400008, 0.00521837, 8.7452e-11, -3.80446e-11, 6.22442e-05, 0.000125236,
      0.000175189, -0.00431697, -0.00231078, 1.53093e-10, -0.000127843, -0.000313879, -0.00431697, -0.00231078,
      1.48378e-10, -0.000269892, -0.000204992, -0.00530849, -0.0252191, 0.0263944, -3.22268e-09, 0.00902234, 0.0117702,
      1.97748e-09, -0.0014616, -0.00190676, -3.46669e-11, -2.56074e-11, -0.000945303, 0.000562213, -0.00133483,
      -1.0781e-11, 1.2152e-09, -0.00167752, 1.01179e-10, -9.08818e-11, 1.32054e-09, 4.95403e-09, -0.00167752,
      -9.16991e-11, 3.90207e-10, -3.87718e-10, -8.02206e-10, -3.22268e-09, 0.00586168, 1.17961e-10, -8.0876e-10,
      -0.0246149, 2.19493e-10, -1.51248e-09, -0.000657354, -0.000503887, -5.7127e-11, -5.47433e-11, -3.87554e-11,
      -0.0190573, -0.00675013, 5.50973e-11, 0.00148885, -0.000705467, 0.00508405, 0.00214883, -3.95918e-11, -0.000871478,
      -0.000718517, 0.000948101, 0.00156873, 0.00902234, 1.17961e-10, 0.0113383, -0.00134147, 8.17922e-11, -0.00942893,
      -0.00424177, 1.0493e-11, 9.31353e-11, -0.000250714, -0.000408121, -0.0001564, 0.000138104, 0.000409342,
      -6.39245e-11, -8.35569e-05, -0.00125455, -0.0183672, -0.00641206, 1.08879e-10, -0.00134043, 0.00110577,
      0.00123686, 0.00204651, 0.0117702, -8.0876e-10, -0.00134147, 0.0106165, -8.21835e-10, -0.00424177, -0.0117111,
      1.29167e-10, -2.16897e-11, 8.89242e-05, -0.000403199, -0.000444212, 3.35572e-10, -6.79241e-10, 0.000483885,
      -1.21659e-10, 8.38774e-11, 1.20741e-09, -2.09349e-09, 0.000483885, 1.60723e-10, -4.27875e-10, -8.56294e-12,
      -8.7714e-10, 1.97748e-09, -0.0246149, 8.17922e-11, -8.21835e-10, 0.0104528, 1.2838e-11, 1.26914e-09, -0.000156385,
      -0.000119876, -5.19244e-11, -1.13761e-11, 3.30078e-11, -0.00623952, 0.00114005, -3.09858e-11, 0.00107889,
      5.94147e-06, 0.00143661, 0.000461806, -2.24067e-11, 0.000217342, 0.000266852, 0.000904666, 0.00400008, -0.0014616,
      2.19493e-10, -0.00942893, -0.00424177, 1.2838e-11, 0.00433404, 0.00117064, 1.29879e-11, -2.81993e-11, -0.00176266,
      -0.00024869, 2.40757e-05, -0.000190833, 0.000784914, 6.48358e-11, 0.000498621, 0.000290727, -0.00607488,
      0.00130481, 9.61837e-11, -0.000225828, 0.00115228, 0.0011802, 0.00521837, -0.00190676, -1.51248e-09, -0.00424177,
      -0.0117111, 1.26914e-09, 0.00117064, 0.00496387, 6.33792e-11, -6.68073e-12, 0.000385743, 0.000509689, -0.00151892,
      -4.00913e-11, 1.25614e-11, -2.41836e-06, 6.05172e-12, 2.69224e-12, 1.21818e-11, -6.3105e-11, 0.000838358,
      -3.90741e-12, -1.45357e-11, 2.3375e-11, 8.7452e-11, -3.46669e-11, -0.000657354, 1.0493e-11, 1.29167e-10,
      -0.000156385, 1.29879e-11, 6.33792e-11, 0.0032637, 1.87206e-05, -5.01137e-11, -4.95967e-12, 1.50133e-11,
      9.93093e-11, -3.2557e-11, 0.000868814, -2.51755e-11, -3.07162e-13, 2.38092e-11, 3.83396e-11, -0.000228034,
      -2.24048e-11, -2.77408e-12, -8.69035e-12, -3.80446e-11, -2.56074e-11, -0.000503887, 9.31353e-11, -2.16897e-11,
      -0.000119876, -2.81993e-11, -6.68073e-12, 1.87206e-05, 0.00325363, -4.49352e-12, 5.5078e-11, -8.48927e-12,
      0.00133348, -0.000219469, -8.87137e-12, -0.000144674, 4.67149e-05, 4.52972e-05, 0.000246225, 2.28275e-12,
      -0.000225506, -9.22014e-05, -5.15423e-05, 6.22442e-05, -0.000945303, -5.7127e-11, -0.000250714, 8.89242e-05,
      -5.19244e-11, -0.00176266, 0.000385743, -5.01137e-11, -4.49352e-12, 0.00280143, 3.78104e-05, 0.000159971,
      -0.000374311, -2.23044e-05, 4.90901e-12, 5.55284e-05, 0.000727652, -0.000774462, 0.000122354, 4.40587e-13,
      0.000592509, 0.000253637, 7.38745e-05, 0.000125236, 0.000562213, -5.47433e-11, -0.000408121, -0.000403199,
      -1.13761e-11, -0.00024869, 0.000509689, -4.95967e-12, 5.5078e-11, 3.78104e-05, 0.0027839, 9.94315e-05,
      0.000513156, 0.000184521, -7.28254e-12, -7.92728e-05, 0.000115256, 0.00125689, -8.43477e-05, 2.03255e-11,
      0.000309832, -6.70315e-05, -4.95276e-05, 0.000175189, -0.00133483, -3.87554e-11, -0.0001564, -0.000444212,
      3.30078e-11, 2.40757e-05, -0.00151892, 1.50133e-11, -8.48927e-12, 0.000159971, 9.94315e-05, 0.00305139;

  EXPECT_NEAR(0.0, (refLocalDensity.array().abs() - localDensity.array().abs()).array().abs().sum(), 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
}

} /* namespace Serenity */

#endif /* POSTHF_MPN_LOCALMP2_TEST_CPP_ */
