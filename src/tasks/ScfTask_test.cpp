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
#include "tasks/ScfTask.h"                       //To be tested.
#include "data/ElectronicStructure.h"            //GetEnergy.
#include "geometry/MolecularSurfaceController.h" //Cavity energy.
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h" //GetElectronicStructure.
#include "testsupply/SystemController__TEST_SUPPLY.h"
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

TEST_F(ScfTaskTest, restricted_cart_ethane) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  Settings settings = sys->getSettings();
  settings.basis.makeSphericalBasis = false;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::HF;
  settings.scf.initialguess = Options::INITIAL_GUESSES::H_CORE;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;

  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, settings);

  ScfTask<RESTRICTED> scfCart(sys);
  scfCart.run();

  double eCart = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  EXPECT_NEAR(-79.1786812939, eCart, 1e-8);

  settings.basis.densFitJ = Options::DENS_FITS::ACD;
  settings.basis.densFitK = Options::DENS_FITS::ACD;
  settings.basis.densFitLRK = Options::DENS_FITS::ACD;
  settings.basis.densFitCorr = Options::DENS_FITS::ACD;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86, settings);

  ScfTask<RESTRICTED> scfChol(sys);
  scfChol.run();

  double eCartChol = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  EXPECT_NEAR(eCart, eCartChol, 1e-4);

  settings.basis.makeSphericalBasis = true;
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP, settings);

  ScfTask<RESTRICTED> scfSph(sys);
  scfSph.run();
  double eSph = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  EXPECT_NEAR(-79.1728833551, eSph, 1e-8);

  EXPECT_NEAR(eCart, eSph, 1e-2);
  SystemController__TEST_SUPPLY::cleanUp();
}

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
  scf.settings.mp2Type = Options::MP2_TYPES::DF;
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
  settings.basis.densFitJ = Options::DENS_FITS::NONE;
  settings.basis.densFitK = Options::DENS_FITS::NONE;
  settings.basis.densFitLRK = Options::DENS_FITS::NONE;
  settings.basis.densFitCorr = Options::DENS_FITS::NONE;
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
  EXPECT_NEAR(-76.4072970518, sys->getElectronicStructure<RESTRICTED>()->getEnergy(), 1e-7);

  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}
TEST_F(ScfTaskTest, Water_PCM_Water) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.cavityFormation = true;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::IEFPCM;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  auto vdwCav = sys->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE_VDW);
  double cavityFormationEnergy = vdwCav->getCavityEnergy();
  EXPECT_NEAR(-76.416394955332379, sys->getElectronicStructure<RESTRICTED>()->getEnergy() - cavityFormationEnergy, 1e-7);
  EXPECT_NEAR(0.007651626114384984, cavityFormationEnergy, 1e-6);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(ScfTaskTest, unrestrictedWater_PCM_Water) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT);
  Settings settings = sys->getSettings();
  settings.pcm.use = true;
  settings.pcm.solvent = Options::PCM_SOLVENTS::WATER;
  settings.pcm.solverType = Options::PCM_SOLVER_TYPES::IEFPCM;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, settings);

  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  EXPECT_NEAR(-76.416394955341488, sys->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-7);
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

TEST_F(ScfTaskTest, unrestricted_DeltaScf) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  ScfTask<UNRESTRICTED> scf(sys);
  scf.run();
  double energy = sys->getElectronicStructure<UNRESTRICTED>()->getEnergy();
  ScfTask<UNRESTRICTED> scf2(sys);
  scf2.settings.exca = {0, 0};
  scf2.run();
  double energy_imom = sys->getElectronicStructure<UNRESTRICTED>()->getEnergy() - energy;
  EXPECT_NEAR(0.289527, energy_imom, 1e-6);

  ScfTask<UNRESTRICTED> scf3(sys);
  scf3.run();
  double energy2 = sys->getElectronicStructure<UNRESTRICTED>()->getEnergy();
  ScfTask<UNRESTRICTED> scf4(sys);
  scf4.settings.excb = {0, 0};
  scf4.settings.momCycles = 100;
  scf4.run();
  double energy_mom = sys->getElectronicStructure<UNRESTRICTED>()->getEnergy() - energy2;
  EXPECT_NEAR(0.289527, energy_mom, 1e-6);

  EXPECT_NEAR(0, energy_imom - energy_mom, 1e-6);
}

TEST_F(ScfTaskTest, restricted_DeltaScf) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs_DFT, true);
  ScfTask<RESTRICTED> scf(sys);
  scf.run();
  double energy = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  ScfTask<RESTRICTED> scf2(sys);
  scf2.settings.exca = {0, 0};
  scf2.run();
  double energy_imom = sys->getElectronicStructure<RESTRICTED>()->getEnergy() - energy;
  EXPECT_NEAR(0.806179, energy_imom, 1e-6);

  ScfTask<RESTRICTED> scf3(sys);
  scf3.run();
  double energy2 = sys->getElectronicStructure<RESTRICTED>()->getEnergy();
  ScfTask<RESTRICTED> scf4(sys);
  scf4.settings.exca = {0, 0};
  scf4.settings.momCycles = 100;
  scf4.run();
  double energy_mom = sys->getElectronicStructure<RESTRICTED>()->getEnergy() - energy2;
  EXPECT_NEAR(0.806179, energy_mom, 1e-6);

  EXPECT_NEAR(0, energy_imom - energy_mom, 1e-6);
}

TEST_F(ScfTaskTest, rohf) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ);
  Settings sysSettings = sys->getSettings();
  sysSettings.scf.rohf = Options::ROHF_TYPES::CUHF;
  sysSettings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, sysSettings);
  sys->setSpin(2);

  // Orca 5.0.3 (June 2023) reference: -113.67554181 Eh

  ScfTask<UNRESTRICTED> cuhf(sys);
  cuhf.run();
  EXPECT_NEAR(-113.6755418076, sys->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-6);

  sysSettings.scf.rohf = Options::ROHF_TYPES::SUHF;
  sysSettings.basis.label = "DEF2-SVP";
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, sysSettings);
  sys->setSpin(2);

  ScfTask<UNRESTRICTED> suhf(sys);
  suhf.run();
  EXPECT_NEAR(-113.6843639808, sys->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-6);
}

TEST_F(ScfTaskTest, CompositeVsCustomFunctional) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, true);

  ScfTask<Options::SCF_MODES::UNRESTRICTED> scf(systemController);
  scf.run();

  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = CompositeFunctionals::XCFUNCTIONALS::B3LYP;
  settings.scfMode = Options::SCF_MODES::UNRESTRICTED;
  settings.basis.label = "6-311G";
#if defined SERENITY_PREFER_XCFUN && defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  settings.customFunc.basicFunctionals = {
      BasicFunctionals::BASIC_FUNCTIONALS::X_SLATER, BasicFunctionals::BASIC_FUNCTIONALS::X_B88_CORR,
      BasicFunctionals::BASIC_FUNCTIONALS::C_LYP, BasicFunctionals::BASIC_FUNCTIONALS::C_VWN};
  settings.customFunc.mixingFactors = {0.80, 0.72, 0.81, 0.19};
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::XCFUN;
#elif defined SERENITY_USE_XCFUN && defined SERENITY_USE_LIBXC
  settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::XC_B3LYP5};
  settings.customFunc.mixingFactors = {1.0};
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::LIBXC;
#elif defined SERENITY_USE_XCFUN
  settings.customFunc.basicFunctionals = {
      BasicFunctionals::BASIC_FUNCTIONALS::X_SLATER, BasicFunctionals::BASIC_FUNCTIONALS::X_B88_CORR,
      BasicFunctionals::BASIC_FUNCTIONALS::C_LYP, BasicFunctionals::BASIC_FUNCTIONALS::C_VWN};
  settings.customFunc.mixingFactors = {0.80, 0.72, 0.81, 0.19};
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::XCFUN;
#else
  settings.customFunc.basicFunctionals = {BasicFunctionals::BASIC_FUNCTIONALS::XC_B3LYP5};
  settings.customFunc.mixingFactors = {1.0};
  settings.customFunc.impl = CompositeFunctionals::IMPLEMENTATIONS::LIBXC;
#endif
  settings.customFunc.hfExchangeRatio = 0.2;
  auto systemController2 =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_311G_3A, settings, 0, 2);

  ScfTask<Options::SCF_MODES::UNRESTRICTED> scf2(systemController2);
  scf2.run();

  EXPECT_NEAR(systemController->getElectronicStructure<UNRESTRICTED>()->getEnergy(),
              systemController2->getElectronicStructure<UNRESTRICTED>()->getEnergy(), 1e-8);
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(systemController2);
}

} /*namespace Serenity*/
