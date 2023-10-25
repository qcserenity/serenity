/**
 * @file LRXPotential_test.cpp
 *
 * @date Jul 12, 2017
 * @author Michael Boeckers
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
#include "potentials/LRXPotential.h"
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

TEST(LRXPotential, CAMB3LYP_restricted) {
  // Perform SCF calculation with CAMB3LYP functional
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  auto es = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  double e0 = es->getEnergyComponentController()->getTotalEnergy();
  EXPECT_NEAR(e0, -76.329704416767839, 1.0e-6);
  auto grad = es->getPotentialBundle()->getGradients();
  // Orca reference values
  EXPECT_NEAR(grad(0, 0), -0.002633230, 1.0e-4);
  EXPECT_NEAR(grad(0, 1), -0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(0, 2), +0.007522045, 1.0e-4);
  EXPECT_NEAR(grad(1, 0), +0.006581398, 1.0e-4);
  EXPECT_NEAR(grad(1, 1), -0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(1, 2), -0.004591275, 1.0e-4);
  EXPECT_NEAR(grad(2, 0), -0.003948168, 1.0e-4);
  EXPECT_NEAR(grad(2, 1), -0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(2, 2), -0.002930770, 1.0e-4);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, CAMB3LYP_unrestricted) {
  // ToDo: Test open-shell system...
  // Perform SCF calculation with CAMB3LYP functional
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP, true);
  double e0 =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergyComponentController()->getTotalEnergy();
  EXPECT_NEAR(e0, -76.329704416738096, 1.0e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, CAMB3LYP_open_shell) {
  Settings settings;
  settings.basis.label = "DEF2-SVP";
  settings.scfMode = UNRESTRICTED;
  settings.name = "TestSystem_WATER_DEF2_SVP";
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP, settings, 0, 2);
  auto es = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  double e0 = es->getEnergyComponentController()->getTotalEnergy();
  // Orca reference value
  EXPECT_NEAR(e0, -76.058271046995, 2.0e-4);
  auto grad = es->getPotentialBundle()->getGradients();
  // Orca reference values (engrad Grid6 NoFinalGrid)
  EXPECT_NEAR(grad(0, 0), +0.000556609, 1.0e-4);
  EXPECT_NEAR(grad(0, 1), +0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(0, 2), -0.101699778, 1.0e-4);
  EXPECT_NEAR(grad(1, 0), -0.098070742, 1.0e-4);
  EXPECT_NEAR(grad(1, 1), +0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(1, 2), +0.026959510, 1.0e-4);
  EXPECT_NEAR(grad(2, 0), +0.097514133, 1.0e-4);
  EXPECT_NEAR(grad(2, 1), +0.000000000, 1.0e-4);
  EXPECT_NEAR(grad(2, 2), +0.074745820, 1.0e-4);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, fock_restricted) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::RESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto lrF = lrxPot.getMatrix();
  EXPECT_NEAR(lrF(0, 0), -4.898955097564e-01, 1.0e-7);
  EXPECT_NEAR(lrF(1, 0), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF(0, 1), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF(1, 1), -4.898955097564e-01, 1.0e-7);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, fock_unrestricted) {
  // ToDo: Test open-shell system...
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::UNRESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto& lrF = lrxPot.getMatrix();
  EXPECT_NEAR(lrF.alpha(0, 0), -4.898955097564e-01, 1.0e-7);
  EXPECT_NEAR(lrF.alpha(1, 0), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF.alpha(0, 1), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF.alpha(1, 1), -4.898955097564e-01, 1.0e-7);
  EXPECT_NEAR(lrF.beta(0, 0), -4.898955097564e-01, 1.0e-7);
  EXPECT_NEAR(lrF.beta(1, 0), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF.beta(0, 1), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF.beta(1, 1), -4.898955097564e-01, 1.0e-7);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, gradients_restricted) {
  // This is only a technical test the actual values are tested in the CAM-B3LYP tests.
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::RESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto grad = lrxPot.getGeomGradients();
  ASSERT_NEAR(grad(0, 0), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(0, 1), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(1, 0), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(1, 1), 0.0, 1.0e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

TEST(LRXPotential, gradients_unrestricted) {
  // This is only a technical test the actual values are tested in the CAM-B3LYP tests.
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::UNRESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto grad = lrxPot.getGeomGradients();
  ASSERT_NEAR(grad(0, 0), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(0, 1), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(1, 0), 0.0, 1.0e-9);
  ASSERT_NEAR(grad(1, 1), 0.0, 1.0e-9);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSystemPath().c_str());
}

} // namespace Serenity
