/**
 * @file ScfTask_test.cpp
 *
 * @author Moritz Bensberg
 * @date Jan 28, 2020
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
#include "tasks/ScfTask.h"                            //To be tested.
#include "data/ElectronicStructure.h"                 //GetEnergy.
#include "settings/Settings.h"                        //Settings.
#include "system/SystemController.h"                  //GetElectronicStructure.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
class ScfTaskTest : public ::testing::Test {
 protected:
  ScfTaskTest() {
  }
  virtual ~ScfTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(ScfTaskTest, restricted_doubleHybridFunctional) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::VERBOSE;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<RESTRICTED> scf(sys);
  scf.run();

  // Changed (12.02.2020) due to separate evaluation of Coulomb and XX
  EXPECT_NEAR(-76.2898331883, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-8);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(ScfTaskTest, unrestricted_doubleHybridFunctional) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<UNRESTRICTED> scf(sys);
  scf.settings.mp2Type = Options::MP2_TYPES::RI;
  scf.run();

  // Changed (12.02.2020) due to separate evaluation of Coulomb and XX
  EXPECT_NEAR(-76.2898331879, sys->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-8);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(ScfTaskTest, restricted_doubleHybridFunctional_LocalMP2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<RESTRICTED> scf(sys);
  scf.settings.mp2Type = Options::MP2_TYPES::LOCAL;
  scf.run();

  // Difference to RIMP2-based result is ~1e-6
  EXPECT_NEAR(-76.289832314507294, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 5e-6);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(ScfTaskTest, restricted_doubleHybridFunctional_Full4CenterMP2) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::MINIMUM;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B2PLYP;
  settings.basis.label = "DEF2-SVP";
  settings.dft.densityFitting = Options::DENS_FITS::NONE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<RESTRICTED> scf(sys);
  scf.settings.mp2Type = Options::MP2_TYPES::AO;
  scf.run();

  EXPECT_NEAR(-76.289764909590303, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-7);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

TEST_F(ScfTaskTest, restartFromConverged) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);

  ScfTask<RESTRICTED> scf(sys);
  scf.settings.restart = true;
  scf.run();
  // Serenity 28.01.2020.
  EXPECT_NEAR(-76.4072970518, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-8);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}
TEST_F(ScfTaskTest, Water_PCM_Water) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  EXPECT_NEAR(-76.416436801427736, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ScfTaskTest, unrestrictedWater_PCM_Water) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  EXPECT_NEAR(-76.416436801427736, sys->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ScfTaskTest, Water_skipSCF) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  ScfTask<RESTRICTED> scf2(sys);
  scf2.settings.skipSCF = true;
  scf2.run();
  EXPECT_NEAR(energy, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-9);
}

TEST_F(ScfTaskTest, Water_ADIIS) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energyDIIS = sys->getElectronicStructure<RESTRICTED>()->getEnergy();

  Settings settings = sys->getSettings();
  settings.scf.initialguess = Options::INITIAL_GUESSES::H_CORE;
  settings.scf.useADIIS = true;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);
  ScfTask<RESTRICTED> scfADIIS(sys);
  scfADIIS.run();
  double energyADIIS = sys->getElectronicStructure<RESTRICTED>()->getEnergy();

  EXPECT_NEAR(energyADIIS, energyDIIS, 1e-7);
  SystemController__TEST_SUPPLY::cleanUp();
}
} /*namespace Serenity*/
