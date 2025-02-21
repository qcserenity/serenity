/**
 * @file ERIPotential_test.cpp
 *
 * @date Dec 1, 2016
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
#include "potentials/ERIPotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/Potential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ERIPotentialTest : public ::testing::Test {
 protected:
  ERIPotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
      systemControllerBP86(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86)) {
  }

  virtual ~ERIPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
  std::shared_ptr<SystemController> systemControllerBP86;
};

/**
 * @test ERIPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */

TEST_F(ERIPotentialTest, H2_FockMatrix) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);
}

/**
 * RESTRICTED FOCK MATRIX
 */

TEST_F(ERIPotentialTest, H2_FockMatrix_CD) {
  Settings settings = systemController->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::CD;
  settings.basis.densFitK = Options::DENS_FITS::CD;
  settings.basis.densFitLRK = Options::DENS_FITS::CD;
  settings.basis.densFitCorr = Options::DENS_FITS::CD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_ACD) {
  Settings settings = systemController->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_ACCD) {
  Settings settings = systemController->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::ACCD;
  settings.basis.densFitK = Options::DENS_FITS::ACCD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_RI) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.584044499, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.038862792, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.534599325, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.038862792, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.534599325, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_RI1) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_RI2) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.581134944, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036392024, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533273657, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036392024, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533273657, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACD1) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACD2) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACCD1) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::ACCD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACCD2) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::ACCD;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

/**
 * UNRESTRICTED FOCK MATRIX
 */

TEST_F(ERIPotentialTest, H2_FockMatrix_CD_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.basis.densFitJ = Options::DENS_FITS::CD;
  settings.basis.densFitK = Options::DENS_FITS::CD;
  settings.basis.densFitLRK = Options::DENS_FITS::CD;
  settings.basis.densFitCorr = Options::DENS_FITS::CD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_ACD_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_ACCD_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.basis.densFitJ = Options::DENS_FITS::ACCD;
  settings.basis.densFitK = Options::DENS_FITS::ACCD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.58365844, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.038637890, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.53449145, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_RI_U) {
  Settings settings = systemController->getSettings();
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::HF;
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.584044499, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.038862792, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.534599325, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.038862792, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.534599325, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.584044499, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.038862792, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.534599325, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.038862792, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.534599325, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_RI1_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::RI;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_RI2_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::RI;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::RI;
  settings.basis.densFitCorr = Options::DENS_FITS::RI;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.581134944, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036392024, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533273657, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036392024, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533273657, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.581134944, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036392024, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533273657, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036392024, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533273657, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACD1_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACD2_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACCD1_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::ACCD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

TEST_F(ERIPotentialTest, H2_FockMatrix_LRX_ACCD2_U) {
  Settings settings = systemController->getSettings();
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::ACCD;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::ACCD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACCD;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings);

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      sys->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
  EXPECT_NEAR(F.alpha(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.alpha(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(F.beta(0, 0), 1.580783084, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 1), 1.036184563, 1e-5);
  EXPECT_NEAR(F.beta(0, 2), 0.533173129, 1e-5);
  EXPECT_NEAR(F.beta(0, 3), 0.0, 1e-5);

  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys);
}

/**
 * GRADIENTS
 */

/**
 * @test ERIPotentialTest
 * @brief Tests the ERI gradient part of an H2
 */
TEST_F(ERIPotentialTest, H2_rgrad) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController, dMat, 1.0, 0.0, 0.0, 0, 0.0);

  auto derivative = coulPot.getGeomGradients();
  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), 0.34677543809435701, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), -0.34677543809435701, 1e-5);
}

/**
 * @test ERIPotentialTest
 * @brief Tests the ERI gradient part of an H2
 */
TEST_F(ERIPotentialTest, H2_ugrad) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemController, dMat, 1.0, 0.0, 0.0, 0.0, 0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), 0.34677543809435701, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), -0.34677543809435701, 1e-5);
}

/**
 * @test ERIPotentialTest
 * @brief Tests the BP86 ERI gradient part of an H2
 */
TEST_F(ERIPotentialTest, H2_grad_rBP86) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat =
      systemControllerBP86->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  Settings settings = systemControllerBP86->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  settings.basis.label = "6-31+G";
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86, settings);

  ERIPotential<Options::SCF_MODES::RESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), 0.71441255550859717, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), -0.71441255550859717, 1e-5);
}

/**
 * @test ERIPotentialTest
 * @brief Tests the BP86 ERI gradient part of an H2
 */
TEST_F(ERIPotentialTest, H2_grad_uBP86) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat =
      systemControllerBP86->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();

  Settings settings = systemControllerBP86->getSettings();
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86, settings);

  ERIPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(sys, dMat, 0.0, 0.0, 0.0, 0.0, 0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), 0.71441255550859717, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), -0.71441255550859717, 1e-5);
}

} // namespace Serenity
