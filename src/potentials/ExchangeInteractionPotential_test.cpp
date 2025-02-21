/**
 * @file ExchangeInteractionPotential_test.cpp
 * @author: Kevin Klahr
 *
 * @date 29. November 2016
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
#include "potentials/ExchangeInteractionPotential.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/CoulombInteractionPotential.h"
#include "potentials/ERIPotential.h"
#include "potentials/Potential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ExchangeInteractionPotentialTest : public ::testing::Test {
 protected:
  ExchangeInteractionPotentialTest()
    : systemControllerAct(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerEnv(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)),
      systemControllerWater(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DISTORTED_MINBAS)) {
    auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, true);
    Settings settings = sys->getSettings();
    settings.basis.densFitJ = Options::DENS_FITS::NONE;
    systemControllerAct =
        SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE, settings);
  }

  virtual ~ExchangeInteractionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerAct;
  std::shared_ptr<SystemController> systemControllerEnv;
  std::shared_ptr<SystemController> systemControllerWater;
};

/**
 * @test HFInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */

TEST_F(ExchangeInteractionPotentialTest, H2_dimer_FockMatrix_Restricted) {
  auto basis = systemControllerAct->getAtomCenteredBasisController();

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemControllerAct, dMat, 1.0,
                                                       basis->getPrescreeningThreshold(), 0.0, 0.0, 0, false);

  FockMatrix<Options::SCF_MODES::RESTRICTED> FNoInt = coulPot.getMatrix();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatAct;
  dMatAct.push_back(systemControllerAct->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ExchangeInteractionPotential<Options::SCF_MODES::RESTRICTED> potential2(basis, dMatAct, 1.0,
                                                                          basis->getPrescreeningThreshold());
  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> potential3(systemControllerAct, {systemControllerEnv},
                                                                         basis, dMatAct);

  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt = potential2.getMatrix();
  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt2 = potential3.getMatrix();
  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt3 = FInt + FInt2;

  EXPECT_NEAR(FNoInt(0, 0), FInt3(0, 0), 1e-5);
  EXPECT_NEAR(FNoInt(1, 0), FInt3(1, 0), 1e-5);
  EXPECT_NEAR(FNoInt(2, 0), FInt3(2, 0), 1e-5);
  EXPECT_NEAR(FNoInt(0, 1), FInt3(0, 1), 1e-5);
  EXPECT_NEAR(FNoInt(0, 2), FInt3(0, 2), 1e-5);
  EXPECT_NEAR(FNoInt(0, 3), FInt3(0, 3), 1e-5);
}

/**
 * @test HFInteractionPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer in the unrestricted case
 */

TEST_F(ExchangeInteractionPotentialTest, H2_dimer_FockMatrix_Unrestricted) {
  auto basis = systemControllerAct->getAtomCenteredBasisController();

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemControllerAct, dMat, 1.0,
                                                         basis->getPrescreeningThreshold(), 0.0, 0.0, 0, false);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FNoInt = coulPot.getMatrix();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatAct;
  dMatAct.push_back(
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());

  ExchangeInteractionPotential<Options::SCF_MODES::UNRESTRICTED> potential2(basis, dMatAct, 1.0,
                                                                            basis->getPrescreeningThreshold());
  CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED> potential3(systemControllerAct, {systemControllerEnv},
                                                                           basis, dMatAct);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt = potential2.getMatrix();
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt2 = potential3.getMatrix();
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt3 = FInt + FInt2;

  for_spin(FNoInt, FInt3) {
    EXPECT_NEAR(FNoInt_spin(0, 0), FInt3_spin(0, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(1, 0), FInt3_spin(1, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(2, 0), FInt3_spin(2, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 1), FInt3_spin(0, 1), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 2), FInt3_spin(0, 2), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 3), FInt3_spin(0, 3), 1e-5);
  };
}

TEST_F(ExchangeInteractionPotentialTest, H2_dimer_FockMatrix_Triplett) {
  systemControllerAct->setSpin(2);
  auto basis = systemControllerAct->getAtomCenteredBasisController();

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemControllerAct, dMat, 1.0,
                                                         basis->getPrescreeningThreshold(), 0.0, 0.0, 0, false);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FNoInt = coulPot.getMatrix();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> dMatAct;
  dMatAct.push_back(
      systemControllerAct->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());

  ExchangeInteractionPotential<Options::SCF_MODES::UNRESTRICTED> potential2(basis, dMatAct, 1.0,
                                                                            basis->getPrescreeningThreshold());
  CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED> potential3(systemControllerAct, {systemControllerEnv},
                                                                           basis, dMatAct);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt = potential2.getMatrix();
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt2 = potential3.getMatrix();
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> FInt3 = FInt + FInt2;

  for_spin(FNoInt, FInt3) {
    EXPECT_NEAR(FNoInt_spin(0, 0), FInt3_spin(0, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(1, 0), FInt3_spin(1, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(2, 0), FInt3_spin(2, 0), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 1), FInt3_spin(0, 1), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 2), FInt3_spin(0, 2), 1e-5);
    EXPECT_NEAR(FNoInt_spin(0, 3), FInt3_spin(0, 3), 1e-5);
  };
}

TEST_F(ExchangeInteractionPotentialTest, H2O_FockMatrix) {
  auto basis = systemControllerWater->getAtomCenteredBasisController();

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerWater->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemControllerWater, dMat, 1.0,
                                                       basis->getPrescreeningThreshold(), 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> FNoInt = coulPot.getMatrix();

  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> dMatAct;
  dMatAct.push_back(
      systemControllerWater->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ExchangeInteractionPotential<Options::SCF_MODES::RESTRICTED> potential2(basis, dMatAct, 1.0,
                                                                          basis->getPrescreeningThreshold());
  CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> potential3(systemControllerAct, {systemControllerEnv},
                                                                         basis, dMatAct);

  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt = potential2.getMatrix();
  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt2 = potential3.getMatrix();
  FockMatrix<Options::SCF_MODES::RESTRICTED> FInt3 = FInt + FInt2;

  EXPECT_NEAR(FNoInt(0, 0), FInt3(0, 0), 1e-5);
  EXPECT_NEAR(FNoInt(1, 0), FInt3(1, 0), 1e-5);
  EXPECT_NEAR(FNoInt(2, 0), FInt3(2, 0), 1e-5);
  EXPECT_NEAR(FNoInt(0, 1), FInt3(0, 1), 1e-5);
  EXPECT_NEAR(FNoInt(0, 2), FInt3(0, 2), 1e-5);
  EXPECT_NEAR(FNoInt(0, 3), FInt3(0, 3), 1e-5);
}

} // namespace Serenity
