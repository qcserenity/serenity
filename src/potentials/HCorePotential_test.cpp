/**
 * @file HCorePotential_test.cpp
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
#include "potentials/HCorePotential.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "geometry/Geometry.h"
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
#include "potentials/Potential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/OrbitalsIOTask.h" // Read MOs from disc.
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class HCorePotentialTest : public ::testing::Test {
 protected:
  HCorePotentialTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
      systemControllerBP86(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86)),
      i2(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE)) {
  }

  virtual ~HCorePotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemController;
  std::shared_ptr<SystemController> systemControllerBP86;
  std::shared_ptr<SystemController> i2;
};

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore part of the Fock Matrix of an H2
 */
TEST_F(HCorePotentialTest, H2_rFockMatrix) {
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(systemController);
  FockMatrix<Options::SCF_MODES::RESTRICTED> F = hCorePot.getMatrix();
  EXPECT_NEAR(F(0, 0), -0.513929914348, 1e-5);
  EXPECT_NEAR(F(1, 0), -1.008815705702, 1e-5);
  EXPECT_NEAR(F(2, 0), -0.669356697423, 1e-5);
  EXPECT_NEAR(F(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(5, 0), -0.1910471699508, 1e-5);
  EXPECT_NEAR(F(6, 0), -0.6329522633973, 1e-5);
  EXPECT_NEAR(F(7, 0), -0.7518036641806, 1e-5);
  EXPECT_NEAR(F(8, 0), -0.5786438117393, 1e-5);
  EXPECT_NEAR(F(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(F(11, 0), 0.5621316727499, 1e-5);
  EXPECT_NEAR(F(1, 1), -1.058271229669, 1e-5);
  EXPECT_NEAR(F(2, 2), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F(3, 3), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F(4, 4), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F(5, 5), 0.1997705812370, 1e-5);
  EXPECT_NEAR(F(6, 6), -0.5139299143488, 1e-5);
  EXPECT_NEAR(F(7, 7), -1.058271229669, 1e-5);
  EXPECT_NEAR(F(8, 8), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F(9, 9), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F(10, 10), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F(11, 11), 0.1997705812370, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore part of the Fock Matrix of an H2
 */
TEST_F(HCorePotentialTest, H2_uFockMatrix) {
  HCorePotential<Options::SCF_MODES::UNRESTRICTED> hCorePot(systemController);
  FockMatrix<Options::SCF_MODES::UNRESTRICTED> F(std::move(hCorePot.getMatrix()));
  EXPECT_NEAR(F.alpha(0, 0), -0.513929914348, 1e-5);
  EXPECT_NEAR(F.alpha(1, 0), -1.008815705702, 1e-5);
  EXPECT_NEAR(F.alpha(2, 0), -0.669356697423, 1e-5);
  EXPECT_NEAR(F.alpha(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(5, 0), -0.1910471699508, 1e-5);
  EXPECT_NEAR(F.alpha(6, 0), -0.6329522633973, 1e-5);
  EXPECT_NEAR(F.alpha(7, 0), -0.7518036641806, 1e-5);
  EXPECT_NEAR(F.alpha(8, 0), -0.5786438117393, 1e-5);
  EXPECT_NEAR(F.alpha(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.alpha(11, 0), 0.5621316727499, 1e-5);
  EXPECT_NEAR(F.alpha(1, 1), -1.058271229669, 1e-5);
  EXPECT_NEAR(F.alpha(2, 2), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F.alpha(3, 3), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.alpha(4, 4), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.alpha(5, 5), 0.1997705812370, 1e-5);
  EXPECT_NEAR(F.alpha(6, 6), -0.5139299143488, 1e-5);
  EXPECT_NEAR(F.alpha(7, 7), -1.058271229669, 1e-5);
  EXPECT_NEAR(F.alpha(8, 8), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F.alpha(9, 9), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.alpha(10, 10), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.alpha(11, 11), 0.1997705812370, 1e-5);
  EXPECT_NEAR(F.beta(0, 0), -0.513929914348, 1e-5);
  EXPECT_NEAR(F.beta(1, 0), -1.008815705702, 1e-5);
  EXPECT_NEAR(F.beta(2, 0), -0.669356697423, 1e-5);
  EXPECT_NEAR(F.beta(3, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(4, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(5, 0), -0.1910471699508, 1e-5);
  EXPECT_NEAR(F.beta(6, 0), -0.6329522633973, 1e-5);
  EXPECT_NEAR(F.beta(7, 0), -0.7518036641806, 1e-5);
  EXPECT_NEAR(F.beta(8, 0), -0.5786438117393, 1e-5);
  EXPECT_NEAR(F.beta(9, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(10, 0), 0.0, 1e-5);
  EXPECT_NEAR(F.beta(11, 0), 0.5621316727499, 1e-5);
  EXPECT_NEAR(F.beta(1, 1), -1.058271229669, 1e-5);
  EXPECT_NEAR(F.beta(2, 2), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F.beta(3, 3), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.beta(4, 4), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.beta(5, 5), 0.1997705812370, 1e-5);
  EXPECT_NEAR(F.beta(6, 6), -0.5139299143488, 1e-5);
  EXPECT_NEAR(F.beta(7, 7), -1.058271229669, 1e-5);
  EXPECT_NEAR(F.beta(8, 8), -0.8086497482821, 1e-5);
  EXPECT_NEAR(F.beta(9, 9), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.beta(10, 10), 0.4455416095329, 1e-5);
  EXPECT_NEAR(F.beta(11, 11), 0.1997705812370, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore gradient part of an H2
 */
TEST_F(HCorePotentialTest, H2_rgrad) {
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(systemController);

  Eigen::MatrixXd derivative = hCorePot.getGeomGradients();
  derivative -= CoreCoreRepulsionDerivative::calculateDerivative(systemController->getGeometry()->getAtoms());

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), -0.861249627698, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), 0.861249627698, 1e-5);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore gradient part of an H2
 */
TEST_F(HCorePotentialTest, H2_ugrad) {
  HCorePotential<Options::SCF_MODES::UNRESTRICTED> hCorePot(systemController);

  Eigen::MatrixXd derivative = hCorePot.getGeomGradients();
  derivative -= CoreCoreRepulsionDerivative::calculateDerivative(systemController->getGeometry()->getAtoms());

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), -0.861249627698, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), 0.861249627698, 1e-5);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the BP86 HCore gradient part of an H2
 */
TEST_F(HCorePotentialTest, H2_grad_rBP86) {
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(systemControllerBP86);

  Eigen::MatrixXd derivative = hCorePot.getGeomGradients();
  derivative -= CoreCoreRepulsionDerivative::calculateDerivative(systemControllerBP86->getGeometry()->getAtoms());

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), -0.963490664845073, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), 0.963490664845073, 1e-5);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the BP86 HCore gradient part of an H2
 */
TEST_F(HCorePotentialTest, H2_grad_uBP86) {
  HCorePotential<Options::SCF_MODES::UNRESTRICTED> hCorePot(systemControllerBP86);

  Eigen::MatrixXd derivative = hCorePot.getGeomGradients();
  derivative -= CoreCoreRepulsionDerivative::calculateDerivative(systemControllerBP86->getGeometry()->getAtoms());

  EXPECT_NEAR(derivative(0, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0, 2), -0.96350268865774491, 1e-5);
  EXPECT_NEAR(derivative(1, 0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1, 2), 0.96350268865774491, 1e-5);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore gradient part of I2 including ECPs.
 */
TEST_F(HCorePotentialTest, I2_ECPs) {
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(i2);

  Eigen::MatrixXd derivative = hCorePot.getGeomGradients();
  derivative -= CoreCoreRepulsionDerivative::calculateDerivative(systemControllerBP86->getGeometry()->getAtoms());

  EXPECT_NEAR(derivative(0, 0), -9.222522228482052, 1e-5);
  EXPECT_NEAR(derivative(0, 1), -22.076897506958474, 1e-5);
  EXPECT_NEAR(derivative(0, 2), -0.510204212847298, 1e-5);
  EXPECT_NEAR(derivative(1, 0), +9.222522228482052, 1e-5);
  EXPECT_NEAR(derivative(1, 1), +22.076897506958474, 1e-5);
  EXPECT_NEAR(derivative(1, 2), +0.510204212847298, 1e-5);
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::I2_Def2_SVP_PBE);
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore energy for a molecule from the WCCR10 set, including ECPs.
 */
TEST_F(HCorePotentialTest, WCCR1010_ECPs) {
  // Please always load. This system is rather large!
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(sys);

  double energy = hCorePot.getEnergy(sys->getElectronicStructure<RESTRICTED>()->getDensityMatrix());

  EXPECT_NEAR(energy, -10893.5666665412, 2e-5); // Compared to turbomole result.
                                                // The accuracy is limited by accuracy of the
  // orbital coefficients on disk.
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore energy for a molecule from the WCCR10 set, including ECPs, and external charges.
 */
TEST_F(HCorePotentialTest, ExternalChargesAndECPs) {
  // Please always load. This system is rather large!
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  std::string pathToTestsResources;
  std::string moPath;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/TestSystem_WCCR1010_P1_Def2_SVP_HF/external_charges.dat";
    moPath = (std::string)env_p + "testresources/TestSystem_WCCR1010_P1_Def2_SVP_HF";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  auto settings = sys->getSettings();
  settings.extCharges.externalChargesFile = pathToTestsResources;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF, settings, 1);
  OrbitalsIOTask<RESTRICTED> read(sys);
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::SERENITY;
  read.settings.path = moPath;
  read.run();

  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(sys);
  double energy = hCorePot.getEnergy(sys->getElectronicStructure<RESTRICTED>()->getDensityMatrix());

  // Taken from "working" implementation (23.01.2025); fixed an error during file reading.
  EXPECT_NEAR(energy, -10893.052566566843, 2e-5);

  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot2(sys);
  double energy2 = hCorePot.getEnergy(sys->getElectronicStructure<RESTRICTED>()->getDensityMatrix());
  EXPECT_NEAR(energy2, energy, 1e-9);

  // orbital coefficients on disk.
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(HCorePotentialTest, Constructor) {
  // Please always load. This system is rather large!
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  std::string pathToTestsResources;
  std::string moPath;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/TestSystem_WCCR1010_P1_Def2_SVP_HF/external_charges.dat";
    moPath = (std::string)env_p + "testresources/TestSystem_WCCR1010_P1_Def2_SVP_HF";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  auto settings = sys->getSettings();
  settings.extCharges.externalChargesFile = pathToTestsResources;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF, settings, 1);

  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(sys);
  EXPECT_NO_THROW(hCorePot.getMatrix());

  // orbital coefficients on disk.
  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test HCorePotentialTest
 * @brief Tests the HCore energy for an H2 molecule, including a constant external potential defined on a grid
 * (RESTRICTED).
 */
TEST_F(HCorePotentialTest, externalGridPotential) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE);
  std::string pathToTestsResources;
  std::string moPath;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/TestSystem_H2_DEF2_TZVP_PBE/externalGridPotential.dat";
    moPath = (std::string)env_p + "testresources/TestSystem_H2_DEF2_TZVP_PBE";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }
  auto settings = sys->getSettings();
  settings.externalGridPotential = pathToTestsResources;
  sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE, settings, 0, 0);
  OrbitalsIOTask<RESTRICTED> read(sys);
  read.settings.fileFormat = Options::ORBITAL_FILE_TYPES::SERENITY;
  read.settings.path = moPath;
  read.run();

  HCorePotential<Options::SCF_MODES::RESTRICTED> hCorePot(sys);
  double energy = hCorePot.getEnergy(sys->getElectronicStructure<RESTRICTED>()->getDensityMatrix());

  // orbital coefficients from disk (calculated without the external potential) lead to this energy. Since the potential
  // used here is 1.0 at each grid point, it should increase the energy by 2 since there are 2 electrons in H2
  EXPECT_NEAR(energy, -2.50308689114 + 2, 1e-6);

  SystemController__TEST_SUPPLY::forget(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_PBE);
  SystemController__TEST_SUPPLY::cleanUp();
}

} // namespace Serenity
