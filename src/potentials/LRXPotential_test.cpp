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
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class LRXPotentialTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(LRXPotentialTest, CAMB3LYP_restricted) {
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
}

TEST_F(LRXPotentialTest, CAMB3LYP_unrestricted) {
  // Perform SCF calculation with CAMB3LYP functional
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP, true);
  double e0 =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergyComponentController()->getTotalEnergy();
  EXPECT_NEAR(e0, -76.329704416738096, 1.0e-6);
}

TEST_F(LRXPotentialTest, CAMB3LYP_open_shell) {
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
}

TEST_F(LRXPotentialTest, fock_restricted) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::RESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto lrF = lrxPot.getMatrix();
  EXPECT_NEAR(lrF(0, 0), -4.898955097564e-01, 1.0e-7);
  EXPECT_NEAR(lrF(1, 0), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF(0, 1), -4.564752841008e-01, 1.0e-7);
  EXPECT_NEAR(lrF(1, 1), -4.898955097564e-01, 1.0e-7);
}

TEST_F(LRXPotentialTest, fock_unrestricted) {
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
}

TEST_F(LRXPotentialTest, fock_unrestricted_openshell) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::UNRESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto& lrF = lrxPot.getMatrix();
  Eigen::MatrixXd refA(6, 6), refB(6, 6);
  refA << -1.096330894, -0.2748540676, 7.25357183089e-18, -0.00275241668942, 5.37362621595e-18, -0.172817650261,
      -0.2748540676, -0.895668655378, 2.27046772122e-17, 0.00366012762939, -6.78779344167e-17, -0.685932870933,
      7.25357183089e-18, 2.27046772122e-17, -0.767221966523, -4.74877307003e-17, 1.1731592037e-16, -3.72849249742e-19,
      -0.00275241668942, 0.00366012762939, -4.74877307003e-17, -0.743174754486, -1.29775089065e-16, -0.347276831139,
      5.37362621595e-18, -6.78779344167e-17, 1.1731592037e-16, -1.29775089065e-16, -0.767221966523, 6.96088155717e-17,
      -0.172817650261, -0.685932870933, -3.72849249742e-19, -0.347276831139, 6.96088155717e-17, -0.704584273113;

  refB << -1.09484827217, -0.265239353662, -1.44397844636e-16, -0.00271521737145, 5.6335119761e-17, -0.166294442379,
      -0.265239353662, -0.81623687606, -3.98829277776e-16, 0.00630742887506, -7.31243021178e-16, -0.633066890845,
      -1.44397844636e-16, -3.98829277776e-16, -0.586412133429, 2.10862436311e-17, 0.289981084824, -4.04722756702e-16,
      -0.00271521737145, 0.00630742887506, 2.10862436311e-17, -0.734052384507, -2.32319550411e-16, -0.342824032005,
      5.6335119761e-17, -7.31243021178e-16, 0.289981084824, -2.32319550411e-16, -0.266256751031, -5.99599509616e-16,
      -0.166294442379, -0.633066890845, -4.04722756702e-16, -0.342824032005, -5.99599509616e-16, -0.669491258416;

  EXPECT_NEAR((lrF.alpha - refA).norm(), 0, 1.0e-7);
  EXPECT_NEAR((lrF.beta - refB).norm(), 0, 1.0e-7);
}

TEST_F(LRXPotentialTest, gradients_restricted) {
  // This is only a technical test the actual values are tested in the CAM-B3LYP tests.
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::RESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto grad = lrxPot.getGeomGradients();
  EXPECT_NEAR(grad(0, 0), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(0, 1), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(1, 0), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(1, 1), 0.0, 1.0e-9);
}

TEST_F(LRXPotentialTest, gradients_unrestricted) {
  // This is only a technical test the actual values are tested in the CAM-B3LYP tests.
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController =
      systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::UNRESTRICTED> lrxPot(systemController, dMatController, 1.0, 1.0e-20, 0.0, 0.0, 0, 1.0);
  auto grad = lrxPot.getGeomGradients();
  EXPECT_NEAR(grad(0, 0), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(0, 1), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(1, 0), 0.0, 1.0e-9);
  EXPECT_NEAR(grad(1, 1), 0.0, 1.0e-9);
}

} // namespace Serenity
