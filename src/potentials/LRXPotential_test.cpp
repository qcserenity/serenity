/**
 * @file LRXPotential_test.cpp
 *
 * @date Jul 12, 2017
 * @author Michael Boeckers
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
#include "data/ElectronicStructure.h"
#include "energies/EnergyComponentController.h"
#include "potentials/LRXPotential.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

TEST(LRXPotential,CAMB3LYP_RESTRICTED) {
  //Perform SCF calculation with CAMB3LYP functional
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  double e0 = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergyComponentController()->getTotalEnergy();
  EXPECT_NEAR(e0,-76.329704416767839,1.0e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  rmdir(systemController->getSettings().path.c_str());
}

TEST(LRXPotential,CAMB3LYP_UNRESTRICTED) {
  //ToDo: Test open-shell system...
  //Perform SCF calculation with CAMB3LYP functional
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  double e0 = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getEnergyComponentController()->getTotalEnergy();
  EXPECT_NEAR(e0,-76.329704416738096,1.0e-6);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WATER_DEF2_SVP_CAMB3LYP);
  rmdir(systemController->getSettings().path.c_str());
}

TEST(LRXPotential,FOCK_RESTRICTED) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::RESTRICTED> lrxPot(dMatController,1.0,1.0e-20,1.0);
  auto lrF=lrxPot.getMatrix();
  EXPECT_NEAR(lrF(0,0),-4.898955097564e-01,1.0e-7);
  EXPECT_NEAR(lrF(1,0),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF(0,1),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF(1,1),-4.898955097564e-01,1.0e-7);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSettings().path.c_str());
}

TEST(LRXPotential,FOCK_UNRESTRICTED) {
  //ToDo: Test open-shell system...
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  auto dMatController = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController();
  LRXPotential<Options::SCF_MODES::UNRESTRICTED> lrxPot(dMatController,1.0,1.0e-20,1.0);
  auto& lrF=lrxPot.getMatrix();
  EXPECT_NEAR(lrF.alpha(0,0),-4.898955097564e-01,1.0e-7);
  EXPECT_NEAR(lrF.alpha(1,0),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF.alpha(0,1),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF.alpha(1,1),-4.898955097564e-01,1.0e-7);
  EXPECT_NEAR(lrF.beta(0,0),-4.898955097564e-01,1.0e-7);
  EXPECT_NEAR(lrF.beta(1,0),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF.beta(0,1),-4.564752841008e-01,1.0e-7);
  EXPECT_NEAR(lrF.beta(1,1),-4.898955097564e-01,1.0e-7);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
  rmdir(systemController->getSettings().path.c_str());
}

}




