/**
 * @file OrbitalController_test.cpp
 *
 * @date Aug 21, 2018
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
#include "data/OrbitalController.h"
#include "data/ElectronicStructure.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class OrbitalControllerTest
 * @brief Sets everything up for the tests of OrbitalController.h/.cpp .
 */
class OrbitalControllerTest : public ::testing::Test {
 protected:
  OrbitalControllerTest() {
  }

  virtual ~OrbitalControllerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests elimination of linear dependencies from transformation matrix.
 */
TEST_F(OrbitalControllerTest, SCF_WithLinearDependendBasisSet) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto H2Linear = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT, true);
  Settings settings = H2Linear->getSettings();
  settings.basis.label = "AUG-CC-PVQZ";
  settings.scf.canOrthThreshold = 5e-4;
  H2Linear = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT, settings);

  ScfTask<SPIN> task(H2Linear);
  task.run();

  auto orbs = H2Linear->getElectronicStructure<SPIN>()->getMolecularOrbitals()->getCoefficients();
  auto eps = H2Linear->getElectronicStructure<SPIN>()->getMolecularOrbitals()->getEigenvalues();

  // Three columns should be removed from the calculation.
  // Accordingly the columns should have been set to zero and the eigenvalues should be inf.
  EXPECT_EQ(orbs.rightCols(3).sum(), 0.0);
  EXPECT_TRUE(eps.tail(3).sum() > 1e+20);
}

/**
 * @test
 * @brief Tests setting and resetting core orbitals (includes ECPs.).
 */
TEST_F(OrbitalControllerTest, setCoreOrbitals) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  auto orbitalController = system->getActiveOrbitalController<RESTRICTED>();
  unsigned int nCore_1 = system->getNCoreElectrons();
  unsigned int nCoreFromOrbitalController = orbitalController->getNCoreOrbitals();
  EXPECT_EQ(nCore_1, nCoreFromOrbitalController * 2);
  orbitalController->setCoreOrbitalsByEnergyCutOff(-5);
  unsigned int nCore_2 = orbitalController->getOrbitalFlags().sum();
  EXPECT_EQ(nCore_1, 82);
  EXPECT_EQ(nCore_2, 41);
  orbitalController->setCoreOrbitalsByEnergyCutOff(-10);
  unsigned int nCore_3 = orbitalController->getOrbitalFlags().sum();
  EXPECT_EQ(nCore_3, 29);
  unsigned int nCoreSet = 40;
  orbitalController->setCoreOrbitalsByNumber(nCoreSet);
  unsigned int nCore_4 = orbitalController->getOrbitalFlags().sum();
  EXPECT_EQ(nCoreSet, nCore_4);
}

/**
 * @test
 * @brief Test core orbitals from initial guess.
 */
TEST_F(OrbitalControllerTest, coreOrbitalsFromGuess) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE, true);
  auto nCoreElectronsInTable = system->getNCoreElectrons();
  auto orbitalController = system->getActiveOrbitalController<RESTRICTED>();
  const unsigned int nCoreOrbitals = orbitalController->getNCoreOrbitals();
  EXPECT_EQ(nCoreElectronsInTable, nCoreOrbitals * 2);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Test get all valence orbitals.
 */
TEST_F(OrbitalControllerTest, getAllValenceOrbitals) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE, true);
  auto orbitalController = system->getActiveOrbitalController<RESTRICTED>();
  orbitalController->setCoreOrbitalsByNumber(system->getNCoreElectrons() / 2);
  orbitalController->setRydbergOrbitalsByNumber(
      system->getBasisController()->getNBasisFunctions() -
      system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION)->getNBasisFunctions());
  auto valenceOrbitals = orbitalController->getAllValenceOrbitalIndices();
  for (const auto& iOrb : valenceOrbitals) {
    EXPECT_GT(iOrb, 7);
    EXPECT_LT(iOrb, 26);
  }
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(system->getSystemPath() + "I_FREE/", "I_FREE");
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
