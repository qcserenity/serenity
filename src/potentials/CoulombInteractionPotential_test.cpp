/**
 * @file CoulombInteractionPotential_test.cpp
 *
 * @date Nov 29, 2016
 * @author: Kevin Klahr
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
#include "potentials/CoulombInteractionPotential.h" //To be tested.
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h" //getDensityMatrixController.
#include "data/matrices/DensityMatrixController.h"
#include "system/SystemController.h"                  //Test systems/getElectronicStructure.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <vector>        //std::vector

namespace Serenity {

class CoulombInteractionPotentialTest : public ::testing::Test {
 protected:
  CoulombInteractionPotentialTest()
    : systemControllerAct(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerEnv(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)) {
  }

  virtual ~CoulombInteractionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerAct;
  std::shared_ptr<SystemController> systemControllerEnv;
};

/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(CoulombInteractionPotentialTest, H2_dimer_RI_rFockMatrix) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(
      1, systemControllerEnv->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));

  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat,
      systemControllerAct->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB), envAuxBas);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulIntPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 0.7382521102598405, 1e-5);
  EXPECT_NEAR(F(1, 0), 0.48242780904563937, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.4443669712229073, 1e-5);
  EXPECT_NEAR(F(0, 1), 0.48242780904563937, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.4443669712229073, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.40826366525686103, 1e-5);
}

/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */
TEST_F(CoulombInteractionPotentialTest, H2_dimer_RI_uFockMatrix) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(
      1, systemControllerEnv->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));

  CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat,
      systemControllerAct->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB), envAuxBas);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(coulIntPot.getMatrix()));

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 0.7284545374785848, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 0.47597394427921969, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.43571056727363722, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 0.47597394427921969, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.43571056727363722, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.40214364052562462, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), 0.7284545374785848, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 0.47597394427921969, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.43571056727363722, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 0.47597394427921969, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.43571056727363722, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.40214364052562462, 1e-5);
}

/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer without RI
 */
TEST_F(CoulombInteractionPotentialTest, H2_dimer_NORI_uFockMatrix) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(1, systemControllerEnv->getBasisController());

  CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat, systemControllerAct->getBasisController(), envAuxBas);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(coulIntPot.getMatrix()));

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 0.72426735790215213, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 0.47520427741987026, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.43005746072535933, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 0.47520427741987026, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.43005746072535933, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.40051779972633761, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), 0.72426735790215213, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 0.47520427741987026, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.43005746072535933, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 0.47520427741987026, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.43005746072535933, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.40051779972633761, 1e-5);
}

/**
 * @test CoulombInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer without RI
 */
TEST_F(CoulombInteractionPotentialTest, H2_dimer_NORI_rFockMatrix) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(1, systemControllerEnv->getBasisController());

  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat, systemControllerAct->getBasisController(), envAuxBas);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulIntPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 0.73364381174774307, 1e-5);
  EXPECT_NEAR(F(1, 0), 0.48157610077843027, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.43804386389948369, 1e-5);
  EXPECT_NEAR(F(0, 1), 0.48157610077843027, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.43804386389948369, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.40645903959382912, 1e-5);
}

TEST_F(CoulombInteractionPotentialTest, H2_dimer_RI_rGradient) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(
      1, systemControllerEnv->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat,
      systemControllerAct->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB), envAuxBas);

  auto ElElInt = coulIntPot.getGeomGradients();

  EXPECT_NEAR(ElElInt(0, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(0, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(0, 2), 0.67468087220634176, 1e-4);
  EXPECT_NEAR(ElElInt(1, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(1, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(1, 2), -0.04300812300513579, 1e-4);
}

TEST_F(CoulombInteractionPotentialTest, H2_dimer_RI_uGradient) {
  auto actBasis = systemControllerAct->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerEnv->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());
  std::vector<std::shared_ptr<BasisController>> envAuxBas(
      1, systemControllerEnv->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB));
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> coulIntPot(
      systemControllerAct, {systemControllerEnv}, actBasis, envDmat,
      systemControllerAct->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB), envAuxBas);

  auto ElElInt = coulIntPot.getGeomGradients();

  EXPECT_NEAR(ElElInt(0, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(0, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(0, 2), 0.67468087220634176, 1e-4);
  EXPECT_NEAR(ElElInt(1, 0), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(1, 1), 0.0, 1e-4);
  EXPECT_NEAR(ElElInt(1, 2), -0.04300812300513579, 1e-4);
}

} // namespace Serenity
