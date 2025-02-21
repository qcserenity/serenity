/**
 * @file DomainOverlapMatrixController_test.cpp
 *
 * @author Moritz Bensberg
 * @date Dec 16, 2019
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
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h" //To be tested.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"            //Access to the overlap integrals
#include "postHF/LocalCorrelation/LocalCorrelationController.h"    //Easy access to local correlation quantities
#include "postHF/MPn/LocalMP2.h"                                   //Run MP2
#include "system/SystemController.h"                               //Operator+
#include "tasks/LocalizationTask.h"                                //Localize
#include "tasks/ScfTask.h"                                         //Run SCF
#include "testsupply/SystemController__TEST_SUPPLY.h"              //Test resources
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class DomainOverlapMatrixControllerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(DomainOverlapMatrixControllerTest, Water) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<scfMode> scfTask(system);
  scfTask.run();
  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.useBPAlgorithm = false;
  lCSettings.pnoThreshold = 3.33e-7;
  lCSettings.pnoCoreScaling = 1.0;
  lCSettings.singlesPNOFactor = 1.0;
  lCSettings.diisStartResidual = 1e+1;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  lCController->initializeSingles();
  LocalMP2 localMp2(lCController);
  localMp2.calculateEnergyCorrection();
  auto orbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  auto overlapController = lCController->getDomainOverlapMatrixController();

  for (const auto& pair : orbitalPairs) {
    // Check identity for pair and singles.
    const Eigen::MatrixXd identityPair = Eigen::MatrixXd::Identity(pair->t_ij.cols(), pair->t_ij.cols());
    double diff = (identityPair - *overlapController->getS(pair, pair)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-10);
    const Eigen::MatrixXd identitySingles_i =
        Eigen::MatrixXd::Identity(pair->singles_i->t_i.rows(), pair->singles_i->t_i.rows());
    diff = (identitySingles_i - *overlapController->getS(pair->singles_i, pair->singles_i)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-10);
    const Eigen::MatrixXd identitySingles_j =
        Eigen::MatrixXd::Identity(pair->singles_j->t_i.rows(), pair->singles_j->t_i.rows());
    diff = (identitySingles_j - *overlapController->getS(pair->singles_j, pair->singles_j)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-10);
    if (pair->i == pair->j) {
      diff = (identityPair - *overlapController->getS(pair, pair->singles_i)).array().abs().sum();
      EXPECT_NEAR(diff, 0.0, 1e-10);
      diff = (identityPair - *overlapController->getS(pair, pair->singles_j)).array().abs().sum();
      EXPECT_NEAR(diff, 0.0, 1e-10);
    }
  }
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(DomainOverlapMatrixControllerTest, Phenol) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OHPhenol_Def2_SVP_Act, true);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::PhenylPhenol_Def2_SVP_Env, true);
  auto system = *act + *env;
  ScfTask<RESTRICTED> scfTask(system);
  scfTask.run();
  LocalizationTask locTask(system);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.settings.splitValenceAndCore = true;
  locTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.pnoSettings = Options::PNO_SETTINGS::LOOSE;
  lCSettings.ccsdPairThreshold = 1e-5;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  lCController->initializeSingles();
  LocalMP2 localMp2(lCController);
  localMp2.calculateEnergyCorrection();
  auto orbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  auto overlapController = lCController->getDomainOverlapMatrixController();
  auto pair = orbitalPairs[0]; // This should be a 1s--1s diagonal pair
  auto singles = pair->singles_i;
  auto singlesK = pair->coupledPairs[1]->getKSingles();
  const Eigen::MatrixXd& s_ij_i = *overlapController->getS(pair, singles);
  Eigen::MatrixXd refS_ij_i(4, 9);
  refS_ij_i << -0.998581, 2.39709e-06, -5.86614e-07, -1.84007e-06, 0.0521186, 0.0109454, 2.0755e-06, -4.15425e-06,
      3.37613e-07, -2.20565e-06, -0.997249, 0.013913, -0.0302773, 5.42894e-06, -1.73242e-07, 0.025286, 0.0494168,
      0.036081, -5.24641e-07, 0.0117534, 0.997938, 0.0271266, -2.70766e-06, -6.1307e-07, 0.0407876, -0.0365267,
      -0.0157472, 1.50723e-06, 0.0261952, 0.0252093, -0.968311, 1.95185e-06, 1.19175e-07, 0.0116162, 0.0900848, -0.229785;
  const Eigen::MatrixXd& s_ij_k = *overlapController->getS(pair, singlesK);
  Eigen::MatrixXd refS_ij_k(4, 25);
  refS_ij_k << -0.235585, -2.95349e-05, -1.69301e-06, -1.35811e-05, -1.43952e-05, 7.04517e-06, -0.163328, -5.21522e-06,
      6.88775e-06, 9.30571e-06, 0.939233, 0.000217183, 9.60251e-05, 0.17301, 1.08332e-05, 0.0130274, 1.50424e-05,
      2.93963e-06, 1.94772e-05, 5.26893e-06, 0.0604506, 0.0113147, -1.03968e-06, -5.68891e-06, -1.33596e-06,
      -9.08501e-06, 0.0653389, 0.145081, -0.231591, 0.0238929, -0.144442, 2.86473e-05, -0.243576, 0.226561, 0.252449,
      -0.00016058, 0.719865, -0.0370637, -2.80195e-05, -0.226059, 3.71746e-05, 0.348279, 0.053626, -0.14642, -0.0439027,
      6.14326e-06, 3.46397e-06, -0.0152693, 0.0622969, -0.0305288, 1.35287e-05, -0.117646, 0.435791, 0.0811803,
      -0.000843235, -0.198768, 1.82576e-05, -0.0335048, 0.0353105, -0.0274156, 3.70135e-05, -0.144096, 0.279427,
      -0.000164374, 0.473015, -5.32939e-05, 0.373825, -0.448943, -0.233274, -0.124237, -3.45161e-06, 1.32408e-06,
      -0.0566214, -0.0335556, 0.0307467, 3.61375e-05, -0.307882, -0.0606691, -0.060099, 0.0776972, -0.0765977,
      2.7373e-05, -0.366944, -0.261238, 0.120228, 4.38019e-05, -0.158757, -0.4318, 0.000198438, 0.045873, -3.86844e-05,
      0.184201, 0.0537572, 0.243374, -0.526036, -1.15913e-06, 2.27649e-06, 0.00848468, 0.0858649, 0.268191;
  EXPECT_EQ(pair->singles_i, pair->singles_j);
  double diff = (s_ij_i.array().abs() - refS_ij_i.array().abs()).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-5);
  diff = (s_ij_k.array().abs() - refS_ij_k.array().abs()).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-5);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
