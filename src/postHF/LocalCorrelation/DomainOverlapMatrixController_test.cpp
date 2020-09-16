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

#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
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

} /* namespace Serenity */
