/**
 * @file QuasiCanonicalPAODomainConstructor_test.cpp
 *
 * @date Jun 26, 2019
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
#include "analysis/PAOSelection/QuasiCanonicalPAODomainConstructor.h" //To be tested.
#include "data/OrbitalController.h"                                   //canonical orbital eigenvalues.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h"   //Test is through the integral transformation.
#include "postHF/LocalCorrelation/LocalCorrelationController.h"       //produce QCPAOConstructor.
#include "system/SystemController.h"                                  //Test systems.
#include "tasks/ScfTask.h"                                            //Clean electronic structure.
#include "testsupply/SystemController__TEST_SUPPLY.h"                 //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class QuasiCanonicalPAODomainConstructorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test QuasiCanonicalPAODomainConstructorTest.Water
 * @brief Tests QC-PAO construction and fock matrix blocks
 *        for a diagonal pair 1|1.
 */
TEST_F(QuasiCanonicalPAODomainConstructorTest, Water) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.singlesPNOFactor = 1.0;
  lCSettings.doiPAOThreshold = 1e-12;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  auto orbitalPairs = localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  localCorrelationController->initializeSingles();
  auto qcPAOConst = localCorrelationController->produceQCPAOConstructor(1.0, 1.0, false);
  Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
      system->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
      localCorrelationController->getMO3CenterIntegralController(), orbitalPairs, qcPAOConst);
  localCorrelationController->initializeSingles();
  auto singles = localCorrelationController->getSingles();
  auto pair = orbitalPairs[2];
  const Eigen::VectorXd canOrbEnergies =
      system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();
  double f_ii = pair->f_ij; // The pair is actually the 11-pair.
  double diff = f_ii - canOrbEnergies(pair->i);
  diff = (pair->f_ab - canOrbEnergies.segment(5, canOrbEnergies.size() - 5)).array().abs().sum();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
  EXPECT_NEAR(diff, 0.0, 1e-6);
}

} /* namespace Serenity */
