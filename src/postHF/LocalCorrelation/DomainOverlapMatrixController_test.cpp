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
  locTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.useBPAlgorithm = false;
  lCSettings.pnoThreshold = 3.33e-7;
  lCSettings.pnoCoreScaling = 1.0;
  lCSettings.diisStartResidual = 1e+1;
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  auto orbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  lCController->buildOrbitalPairCouplingMap();
  lCController->initializeSingles();
  auto overlapController = lCController->getDomainOverlapMatrixController();

  for (const auto& pair : orbitalPairs) {
    // Check identity for pair and singles.
    const Eigen::MatrixXd identityPair = Eigen::MatrixXd::Identity(pair->t_ij.cols(), pair->t_ij.cols());
    double diff = (identityPair - *overlapController->getS(pair, pair)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-12);
    const Eigen::MatrixXd identitySingles_i =
        Eigen::MatrixXd::Identity(pair->singles_i->t_i.rows(), pair->singles_i->t_i.rows());
    diff = (identitySingles_i - *overlapController->getS(pair->singles_i, pair->singles_i)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-12);
    const Eigen::MatrixXd identitySingles_j =
        Eigen::MatrixXd::Identity(pair->singles_j->t_i.rows(), pair->singles_j->t_i.rows());
    diff = (identitySingles_j - *overlapController->getS(pair->singles_j, pair->singles_j)).array().abs().sum();
    EXPECT_NEAR(diff, 0.0, 1e-12);
    if (pair->i == pair->j) {
      diff = (identityPair - *overlapController->getS(pair, pair->singles_i)).array().abs().sum();
      EXPECT_NEAR(diff, 0.0, 1e-12);
      diff = (identityPair - *overlapController->getS(pair, pair->singles_j)).array().abs().sum();
      EXPECT_NEAR(diff, 0.0, 1e-12);
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
  auto lCController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  auto orbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  lCController->buildOrbitalPairCouplingMap();
  lCController->initializeSingles();
  LocalMP2 localMp2(lCController);
  localMp2.calculateEnergyCorrection();
  auto overlapController = lCController->getDomainOverlapMatrixController();

  auto pair = orbitalPairs[0]; // This should be a 1s--1s diagonal pair
  auto singles = pair->singles_i;
  auto kSet = pair->coupledPairs[1];
  auto singlesK = kSet->getKSingles();
  auto ikPair = kSet->getIKPair();
  auto kjPair = kSet->getKJPair();
  const Eigen::MatrixXd& s_ij_i = *overlapController->getS(pair, singles);
  Eigen::MatrixXd refS_ij_i(4, 9);
  refS_ij_i << -0.998581, 2.39709e-06, -5.86614e-07, -1.84007e-06, 0.0521186, 0.0109454, 2.0755e-06, -4.15425e-06,
      3.37613e-07, -2.20565e-06, -0.997249, 0.013913, -0.0302773, 5.42894e-06, -1.73242e-07, 0.025286, 0.0494168,
      0.036081, -5.24641e-07, 0.0117534, 0.997938, 0.0271266, -2.70766e-06, -6.1307e-07, 0.0407876, -0.0365267,
      -0.0157472, 1.50723e-06, 0.0261952, 0.0252093, -0.968311, 1.95185e-06, 1.19175e-07, 0.0116162, 0.0900848, -0.229785;
  const Eigen::MatrixXd& s_ij_k = *overlapController->getS(pair, singlesK);
  Eigen::MatrixXd refS_ij_k(4, 9);
  refS_ij_k << 0.0757772, -9.16897e-06, -4.41625e-06, 4.50696e-06, -0.0481341, -0.0382635, -5.11354e-06, -2.60864e-06,
      -8.13518e-06, 1.23626e-06, 0.162048, 0.0744579, 0.027375, -1.5673e-06, -1.009e-06, 0.0771397, -0.069708, 0.114236,
      4.03398e-07, 0.127736, 0.112948, -0.0225495, -6.77293e-07, -4.19632e-07, 0.0577649, -0.0170364, 0.161787,
      -1.93608e-06, -0.0318184, -0.0110211, 0.0362626, 1.62839e-06, 1.3951e-06, -0.0138243, -0.0322624, -0.0315281;
  const Eigen::MatrixXd& s_i_k = *overlapController->getS(singles, singlesK);
  Eigen::MatrixXd refS_i_k(9, 9);
  refS_i_k << -0.0782804, 8.11572e-06, 3.81651e-06, -4.46644e-06, 0.0531892, 0.0417046, 4.41417e-06, 3.06305e-06,
      7.25545e-06, -1.1772e-06, -0.153006, -0.0682179, -0.0290041, 1.70775e-06, 1.13192e-06, -0.071222, 0.0695463,
      -0.102282, 1.885e-07, 0.138124, 0.119285, -0.021046, -3.06809e-07, -2.13536e-07, 0.0649448, -0.0230992, 0.170654,
      2.03335e-06, 0.0244322, 0.00697683, -0.0440853, -1.59988e-06, -2.35536e-06, 0.0120508, 0.0485664, 0.0244871,
      -0.040438, 2.02379e-06, 7.82414e-07, -9.13903e-07, 0.0882196, 0.0533854, 1.25909e-07, 3.92624e-06, 2.95397e-06,
      -0.0259928, 5.5693e-08, 5.3068e-07, -2.28608e-06, 0.0348798, 0.0547815, 3.36286e-06, 4.57105e-06, -1.06507e-06,
      2.18115e-06, -0.0633543, -0.0401532, -0.0333424, -4.20777e-06, -1.5107e-06, -0.0491684, 0.091864, -0.0206747,
      2.28345e-06, 0.13843, 0.0819117, -0.0392334, -2.88826e-06, -4.18576e-06, 0.103919, 0.016869, 0.152313,
      -1.21394e-07, 0.0842913, 0.0539549, 0.0052823, -4.12362e-07, 2.99783e-06, 0.0466404, -0.0476044, 0.099748;
  EXPECT_EQ(pair->singles_i, pair->singles_j);
  double diff = (s_ij_i.array().abs() - refS_ij_i.array().abs()).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-5);
  diff = (s_ij_k.array().abs() - refS_ij_k.array().abs()).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-5);
  diff = (s_i_k.array().abs() - refS_i_k.array().abs()).array().abs().sum();
  EXPECT_NEAR(diff, 0.0, 1e-5);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
