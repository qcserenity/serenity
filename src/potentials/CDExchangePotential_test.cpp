/**
 * @file CDExchangePotential_test.cpp
 *
 * @date Apr 04, 2020
 * @author: Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "potentials/CDExchangePotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "integrals/decomposer/TwoElecFourCenterIntDecomposer.h"
#include "potentials/ExchangePotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class CDExchangePotentialTest : public ::testing::Test {
 protected:
  virtual ~CDExchangePotentialTest() = default;

  virtual void TearDown() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test CDCoulombPotentialTest
 * @brief Tests the CHolesky Coulomb Fock Matrix of H2
 */
TEST_F(CDExchangePotentialTest, H2_FockMatrix) {
  auto sysCont = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  Settings settings = sysCont->getSettings();
  settings.basis.densityFitting = Options::DENS_FITS::ACD;
  settings.basis.cdThreshold = 1e-6;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings, 0, 0);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  auto dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()->getDensityMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fxRef(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fxRef->setZero();

  ExchangePotential<Options::SCF_MODES::RESTRICTED> xPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0,
      0.0, 0, 0, 0);

  *fxRef = xPot.getMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fx(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fx->setZero();

  CDExchangePotential<Options::SCF_MODES::RESTRICTED> cdXPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0, 0.0);

  *fx = cdXPot.getMatrix();

  Eigen::MatrixXd diff = (*fx) - (*fxRef);

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-4);

  systemController->getCDIntegralController()->cleanup();
}

/**
 * @test CDCoulombPotentialTest
 * @brief Tests the CHolesky Coulomb Fock Matrix of H2
 */
TEST_F(CDExchangePotentialTest, H2_FockMatrix_accd) {
  auto sysCont = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  Settings settings = sysCont->getSettings();
  settings.basis.densityFitting = Options::DENS_FITS::ACCD;
  settings.basis.cdThreshold = 1e-6;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings, 0, 0);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  auto dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()->getDensityMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fxRef(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fxRef->setZero();

  ExchangePotential<Options::SCF_MODES::RESTRICTED> xPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0,
      0.0, 0, 0, 0);

  *fxRef = xPot.getMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fx(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fx->setZero();

  CDExchangePotential<Options::SCF_MODES::RESTRICTED> cdXPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0, 0.0);

  *fx = cdXPot.getMatrix();

  Eigen::MatrixXd diff = (*fx) - (*fxRef);

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-4);

  systemController->getCDIntegralController()->cleanup();
}

/**
 * @test CDCoulombPotentialTest
 * @brief Tests the CHolesky Coulomb Fock Matrix of H2
 */
TEST_F(CDExchangePotentialTest, F_FockMatrix) {
  auto sysCont = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs);
  Settings settings = sysCont->getSettings();
  settings.basis.densityFitting = Options::DENS_FITS::ACD;
  settings.basis.cdThreshold = 1e-6;
  settings.basis.label = "DEF2-TZVP";
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs, settings, 1, 0);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  auto dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()->getDensityMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fxRef(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fxRef->setZero();

  ExchangePotential<Options::SCF_MODES::RESTRICTED> xPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0,
      0.0, 0, 0, 0);

  *fxRef = xPot.getMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fx(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fx->setZero();

  CDExchangePotential<Options::SCF_MODES::RESTRICTED> cdXPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0, 0.0);

  *fx = cdXPot.getMatrix();

  Eigen::MatrixXd diff = (*fx) - (*fxRef);

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-4);

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  systemController->getCDIntegralController()->cleanup();
}

/**
 * @test CDCoulombPotentialTest
 * @brief Tests the CHolesky Coulomb Fock Matrix of H2
 */
TEST_F(CDExchangePotentialTest, F_FockMatrix_spherical) {
  auto sysCont = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs);
  Settings settings = sysCont->getSettings();
  settings.basis.densityFitting = Options::DENS_FITS::ACD;
  settings.basis.cdThreshold = 1e-6;
  settings.basis.label = "DEF2-TZVP";
  settings.basis.makeSphericalBasis = true;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs, settings, 1, 0);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();

  auto dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController()->getDensityMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fxRef(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fxRef->setZero();

  ExchangePotential<Options::SCF_MODES::RESTRICTED> xPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0,
      0.0, 0, 0, 0);

  *fxRef = xPot.getMatrix();

  FockMatrix<Options::SCF_MODES::RESTRICTED>* fx(
      new FockMatrix<Options::SCF_MODES::RESTRICTED>(systemController->getBasisController()));
  fx->setZero();

  CDExchangePotential<Options::SCF_MODES::RESTRICTED> cdXPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController(), 1.0, 0.0);

  *fx = cdXPot.getMatrix();

  Eigen::MatrixXd diff = (*fx) - (*fxRef);

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-4);

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  systemController->getCDIntegralController()->cleanup();
}

/**
 * @test CDCoulombPotentialTest
 * @brief Tests the CHolesky Coulomb Fock Matrix of H2
 */
TEST_F(CDExchangePotentialTest, H2_FockMatrix_unrestricted) {
  auto sysCont = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  Settings settings = sysCont->getSettings();
  settings.basis.densityFitting = Options::DENS_FITS::ACD;
  settings.basis.cdThreshold = 1e-6;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP, settings, 0, 0);
  systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();

  auto dMat =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController()->getDensityMatrix();

  FockMatrix<Options::SCF_MODES::UNRESTRICTED>* fxRef(
      new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(systemController->getBasisController()));
  fxRef->alpha.setZero();
  fxRef->beta.setZero();

  ExchangePotential<Options::SCF_MODES::UNRESTRICTED> xPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController(), 1.0,
      0.0, 0, 0, 0);

  *fxRef = xPot.getMatrix();

  FockMatrix<Options::SCF_MODES::UNRESTRICTED>* fx(
      new FockMatrix<Options::SCF_MODES::UNRESTRICTED>(systemController->getBasisController()));
  fx->alpha.setZero();
  fx->beta.setZero();

  CDExchangePotential<Options::SCF_MODES::UNRESTRICTED> cdXPot(
      systemController,
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController(), 1.0, 0.0);

  *fx = cdXPot.getMatrix();

  Eigen::MatrixXd diffa = fx->alpha - fxRef->alpha;
  EXPECT_NEAR(diffa.maxCoeff(), 0.0, 1e-4);

  Eigen::MatrixXd diffb = fx->beta - fxRef->beta;
  EXPECT_NEAR(diffb.maxCoeff(), 0.0, 1e-4);

  systemController->getCDIntegralController()->cleanup();
}

} // namespace Serenity
