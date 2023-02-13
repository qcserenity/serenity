/**
 * @file PNOConstructor_test.cpp
 *
 * @date Apr 12, 2019
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

#ifndef POSTHF_PNOCONSTRUCTOR_TEST_CPP_
#define POSTHF_PNOCONSTRUCTOR_TEST_CPP_

/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PNOConstructor.h"     //This will be tested.
#include "data/OrbitalController.h"                   //References and orbital localization.
#include "data/OrbitalPair.h"                         //Orbital pair definition.
#include "data/matrices/FockMatrix.h"                 //Fock matrix definition.
#include "postHF/MPn/LocalMP2.h"                      //Lazy way of generating PNOs.
#include "system/SystemController.h"                  //System controller definition.
#include "tasks/LocalizationTask.h"                   //Orbital localization.
#include "tasks/ScfTask.h"                            //Clean electronic structure from fresh SCF.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test resources.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class PNOConstructorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test PNOConstructorTest.WaterDimer
 * @brief Tests PNO truncation by checking the surviving number of PNOs.
 */
TEST_F(PNOConstructorTest, WaterDimer) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  auto waterDimer = *act1 + *act2;
  ScfTask<scfMode> scfTask(waterDimer);
  scfTask.run();
  LocalizationTask locTask(waterDimer);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.useBPAlgorithm = false;
  lCSettings.pnoThreshold = 1e-8;
  lCSettings.pnoCoreScaling = 1.0;
  lCSettings.diisStartResidual = 1e+1;
  auto lCController = std::make_shared<LocalCorrelationController>(waterDimer, lCSettings);
  const Eigen::MatrixXd f = lCController->getFockMatrix();
  LocalMP2 localMP2(lCController);
  double lMP2Energy = localMP2.calculateEnergyCorrection().sum();
  auto orbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  (void)lMP2Energy;
  Eigen::VectorXi nPNOs = Eigen::VectorXi::Zero(16);
  nPNOs << 7, 4, 17, 20, 11, 20, 17, 21, 20, 13, 20, 19, 16, 21, 21, 20;

  Eigen::MatrixXd c = waterDimer->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  Eigen::MatrixXd f_MO = c.transpose() * f * c;

  unsigned int iPair = 0;
  for (const auto& pair : orbitalPairs) {
    double diff = pair->f_ij - f_MO(pair->i, pair->j);
    EXPECT_NEAR(diff, 0.0, 1e-12);
    if (pair->type != OrbitalPairTypes::VERY_DISTANT) {
      EXPECT_EQ(pair->t_ij.cols(), nPNOs[iPair]);
    }
    else {
      std::cout << "test" << std::endl;
      EXPECT_EQ(0, nPNOs[iPair]);
    }
    iPair++;
    if (iPair >= nPNOs.size())
      break;
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(waterDimer);
  SystemController__TEST_SUPPLY::cleanUp();
}
} /* namespace Serenity */

#endif /* POSTHF_MPN_LOCALMP2_TEST_CPP_ */
