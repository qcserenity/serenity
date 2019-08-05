/**
 * @file Scf_test.cpp
 * 
 * @date Aug 28, 2015
 * @author Thomas Dresselhaus
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
#include "tasks/ScfTask.h"
#include "settings/Settings.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {
/*
 * TODO the resulting numbers are simply compared to old output. Some more proper reference should
 * be used
 */
//===================================
//          RHF
//===================================
/**
 * @test
 * @brief Restricted Hartree--Fock of H2 in a minimal basis using default settings
 */
TEST(BasicSCFTests, H2_MinimalBasis_RHF) {
  Settings settings;
  settings.basis.label = "STO-6G";
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(result->getEnergy(), -1.1253243671827655, 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
/**
 * @test
 * @brief Restricted Hartree--Fock of H2 in a minimal basis using ADIIS.
 */
TEST(BasicSCFTests, H2_MinimalBasis_RHF_ADIIS) {
  Settings settings;
  settings.basis.label = "STO-6G";
  settings.scf.useADIIS = true;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(result->getEnergy(), -1.1253243671827655, 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}

/**
 * @test
 * @brief Restricted Hartree--Fock of H2 in a minimal basis using ADIIS.
 */
TEST(BasicSCFTests, H2O_WaterMonOne_MINBASIS_RHF_ADIIS) {
  Settings settings;
  settings.basis.label = "STO-6G";
  settings.scf.useADIIS = true;
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(result->getEnergy(), -75.679784634374329, 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
//===================================
//          UHF
//===================================
/**
 * @test
 * @brief Unrestricted Hartree--Fock of H2 in a minimal basis using default settings
 */
TEST(BasicSCFTests, H2_MinimalBasis_UHF) {
  Settings settings;
  settings.basis.label = "STO-6G";
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(result->getEnergy(), -1.1253243671827655, 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}

//===================================
//          RKS
//===================================
/**
 * @test 
 * @brief H2 in a minimal basis, DFT, LDA functional (restricted) with RI
 */
TEST(BasicSCFTests, H2_MinimalBasis_RKS_LDA_RI) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = Options::XCFUNCTIONALS::LDA;
  settings.dft.densityFitting = Options::DENS_FITS::RI;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(result->getEnergy(), -1.1299388663978864, 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
/**
 * @test
 * @brief H2 in a minimal basis, DFT, LDA functional (restricted) without RI
 */
TEST(BasicSCFTests, H2_MinimalBasis_RKS_LDA_noRI) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = Options::XCFUNCTIONALS::LDA;
  settings.dft.densityFitting = Options::DENS_FITS::NONE;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(-1.1298876147476988, result->getEnergy(), 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.res.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
//===================================
//          UKS
//===================================
/**
 * @test
 * @brief H2 in a minimal basis, DFT, LDA functional (unrestricted) with RI
 */
TEST(BasicSCFTests, H2_MinimalBasis_UKS_LDA_RI) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = Options::XCFUNCTIONALS::LDA;
  settings.dft.densityFitting = Options::DENS_FITS::RI;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(-1.1299388666806995,result->getEnergy(), 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
/**
 * @test
 * @brief H2 in a minimal basis, DFT, LDA functional (unrestricted) without RI
 */
TEST(BasicSCFTests, H2_MinimalBasis_UKS_LDA_noRI) {
  Settings settings;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = Options::XCFUNCTIONALS::LDA;
  settings.dft.densityFitting = Options::DENS_FITS::NONE;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(-1.1298876147476988, result->getEnergy(), 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
/**
 * @test 
 * @brief H2 triplet in a minimal basis, DFT, LDA functional (unrestricted)
 */
TEST(BasicSCFTests, H2_MinimalBasis_LDA_Triplet) {
  Settings settings;
  settings.spin=2;
  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
  settings.dft.functional = Options::XCFUNCTIONALS::LDA;
  settings.dft.densityFitting = Options::DENS_FITS::NONE;
  settings.basis.label = "STO-6G";
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
  auto result = systemController->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto path=systemController->getSettings().path;
  EXPECT_NEAR(-1.1298876147476988, result->getEnergy(), 1e-6);
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.energies.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.dmat.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.orbs.unres.h5").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.settings").c_str()));
  EXPECT_EQ(0,std::remove((path+"SomeTestSystem.xyz").c_str()));
  std::remove((path).c_str());
}
/**
 * @test
 * @brief H2 in minimal basis, DFT, BLYP functional (restricted)
 */
//TEST(BasicSCFTests, H2_MinimalBasis_BLYP) {
//  Settings settings;
//  settings.method = Options::ELECTRONIC_STRUCTURE_THEORIES::DFT;
//  settings.dft.functional = Options::XCFUNCTIONALS::BLYP;
//  settings.basisLabel[Options::BASIS_PURPOSES::DEFAULT] = "STO-6G";
//  auto systemController = SystemController__TEST_SUPPLY::getSystemController(
//      TEST_SYSTEM_CONTROLLERS::H2_MINBAS, settings);
//  auto result = systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
//  EXPECT_NEAR(result->getEnergy(), -1.1640183820 , 1e-6); // ORCA -1.154908556998
//  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_MINBAS);
//}


// TODO test something with p- and d-functions

} /* namespace Serenity */
