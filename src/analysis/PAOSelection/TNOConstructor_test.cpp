/**
 * @file TNOConstructor_test.cpp
 *
 * @author Moritz Bensberg
 * @date Nov 6, 2019
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
#include "analysis/PAOSelection/TNOConstructor.h"               //To be tested.
#include "data/ElectronicStructure.h"                           //Get Fock matrix.
#include "data/OrbitalController.h"                             //Compare to canonical eigenvalues.
#include "data/OrbitalTriple.h"                                 //Triples definition.
#include "data/PAOController.h"                                 //Transformation of the Fock matrix to TNO basis.
#include "data/matrices/FockMatrix.h"                           //Fock matrix definition.
#include "io/FormattedOutputStream.h"                           //Print level.
#include "postHF/CC/DLPNO_CCSD.h"                               //Initial doubles amplitudes.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //produceTNOConstructor.
#include "system/SystemController.h"                            //Test systems.
#include "tasks/ScfTask.h"                                      //Clean electronic structure.
#include "testsupply/SystemController__TEST_SUPPLY.h"           //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
namespace Serenity {
class TNOConstructorTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test TNOConstructorTest.Water
 * @brief Tests TNO truncation by checking the surviving number of PNOs.
 */
TEST_F(TNOConstructorTest, Water) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP, true);
  ScfTask<scfMode> scfTask(system);
  scfTask.run();

  LocalCorrelationSettings lCSettings;
  lCSettings.useBPAlgorithm = false;
  // Complete PNO and TNO domains.
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.tnoThreshold = 1e-14;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-4, 100);
  localCCSD.calculateElectronicEnergyCorrections().sum();

  auto tnoConstructor = localCorrelationController->produceTNOConstructor();
  auto triples = localCorrelationController->getOrbitalTriples();
  for (auto triple : triples) {
    tnoConstructor->transformToTNOBasis(triple);
  }                         // for triple
  auto triple = triples[0]; // ijk=100, ij=10, ik=10, jk=00
  // Test redundancies
  double diff = (triple->getS_ij_ijk() - triple->getS_ik_ijk()).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-12);
  diff = (triple->getS_ij_ijk() - triple->getS_il_ijk()[0]).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-12);
  // These should be complete TNO domains. Thus, the TNO-pseudo eigenvalues should be
  // identical to the canonical ones.
  auto eigenvalues = system->getActiveOrbitalController<scfMode>()->getEigenvalues();
  diff = (triple->getTNOEigenvalues() - eigenvalues.segment(5, 19)).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-6);

  // Check TNO transformation. The TNOs have to diagonalize the Fock matrix and give
  // the TNO pseudo eigenvalues on the diagonal.
  auto f = system->getElectronicStructure<RESTRICTED>()->getFockMatrix();
  Eigen::MatrixXd p_ijk = localCorrelationController->getPAOController()->getPAOsFromDomain(triple->getPAODomain()) *
                          triple->getTNOCoefficients();
  Eigen::VectorXd f_TNO_diag = (p_ijk.transpose() * f * p_ijk).diagonal();
  diff = (triple->getTNOEigenvalues() - f_TNO_diag).array().abs().sum();
  EXPECT_NEAR(0.0, diff, 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

} /* namespace Serenity */
