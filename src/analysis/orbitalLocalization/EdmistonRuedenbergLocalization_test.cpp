/**
 * @file EdmistonRuedenbergLocalization_test.cpp
 *
 * @date Nov 3, 2016
 * @author David Schnieders
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
#include "analysis/orbitalLocalization/EdmistonRuedenbergLocalization.h"
#include "data/OrbitalController.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/Matrix.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class EdmistonRuedenbergTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test ER Localization Restricted CO minimal basis
 * TODO check whether the density matrix stays the same
 */
TEST_F(EdmistonRuedenbergTest, testLocalizationRestricted) {
  // Create the test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  // Perform SCF
  system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  // Create a copy of the SCF orbitals
  auto orbitals = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(
      *system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>());
  auto coefficients = orbitals->getCoefficients();

  const auto& oneIntController = system->getOneElectronIntegralController();
  const auto& overlaps = oneIntController->getOverlapIntegrals();
  /*
   * Calculate a simple Pipek-Mezey-like localization measure to check whether the algorithm does
   * localize something.
   */
  double convMeasure = 0.0;
  for (unsigned int i = 0; i < 1; ++i) {
    double sumCoeffsThisOrbFirstAtom = 0.0;
    for (unsigned int mu = 0; mu < 2; ++mu) {
      for (unsigned int nu = 0; nu < 4; ++nu) {
        sumCoeffsThisOrbFirstAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasure += sumCoeffsThisOrbFirstAtom * sumCoeffsThisOrbFirstAtom;
    double sumCoeffsThisOrbSecondAtom = 0.0;
    for (unsigned int mu = 2; mu < 4; ++mu) {
      for (unsigned int nu = 0; nu < 4; ++nu) {
        sumCoeffsThisOrbSecondAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasure += sumCoeffsThisOrbSecondAtom * sumCoeffsThisOrbSecondAtom;
  }
  EdmistonRuedenbergLocalization<Options::SCF_MODES::RESTRICTED> localizer(system);
  /*
   * Perform localization
   */
  std::vector<unsigned int> range = {0, 1, 2, 3, 4, 5, 6};
  localizer.localizeOrbitals(*orbitals, 25, range);
  /*
   * check if MOs still orthogonal
   */
  for_spin(coefficients) {
    Matrix<double> overlapMatrix = coefficients_spin.transpose() * overlaps * coefficients_spin;
    Matrix<double> identity(orbitals->getNOrbitals(), orbitals->getNOrbitals());
    identity.setIdentity();
    Matrix<double> diffMatrix = overlapMatrix - identity;
    double diff = fabs(diffMatrix.maxCoeff());
    diff = fabs(diffMatrix.minCoeff()) > diff ? fabs(diffMatrix.minCoeff()) : diff;
    EXPECT_NEAR(0.0, diff, 2.0e-9);
  };
  /*
   * Calculate the convergence measure again for the now localized, occupied orbitals
   * and compare with the value before.
   */
  double convMeasureLoc = 0.0;
  for (unsigned int i = 0; i < 7; ++i) {
    double sumCoeffsThisOrbFirstAtom = 0.0;
    for (unsigned int mu = 0; mu < 5; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        sumCoeffsThisOrbFirstAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasureLoc += sumCoeffsThisOrbFirstAtom * sumCoeffsThisOrbFirstAtom;
    double sumCoeffsThisOrbSecondAtom = 0.0;
    for (unsigned int mu = 5; mu < 10; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        sumCoeffsThisOrbSecondAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasureLoc += sumCoeffsThisOrbSecondAtom * sumCoeffsThisOrbSecondAtom;
  }
  /*
   * TODO the raise in the localization measure is actually much higher. However, due to numerical
   * instabilities in parallel runs the correct result is often not reached.
   */
  EXPECT_GT(convMeasureLoc, convMeasure + 0.1);
}

// Test ER Localization Unrestricted
TEST_F(EdmistonRuedenbergTest, testLocalizationUnrestricted) {
  // Create the test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  // Perform SCF
  system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  // Create a copy of the SCF orbitals
  auto orbitals = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(
      *system->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>());
  auto coefficients = orbitals->getCoefficients();

  EdmistonRuedenbergLocalization<Options::SCF_MODES::UNRESTRICTED> localizer(system);
  /*
   * Perform localization
   */
  auto nOcc = system->getNOccupiedOrbitals<RESTRICTED>();
  std::vector<unsigned int> range;
  for (unsigned int i = 0; i < nOcc; ++i)
    range.push_back(i);
  localizer.localizeOrbitals(*orbitals, 25, range);
  /*
   * check if MOs still orthogonal
   */
  const auto& oneIntController = system->getOneElectronIntegralController();
  const auto& overlaps = oneIntController->getOverlapIntegrals();
  for_spin(coefficients) {
    Matrix<double> overlapMatrix = coefficients_spin.transpose() * overlaps * coefficients_spin;
    Matrix<double> identity(orbitals->getNOrbitals(), orbitals->getNOrbitals());
    identity.setIdentity();
    Matrix<double> diffMatrix = overlapMatrix - identity;
    double diff = fabs(diffMatrix.maxCoeff());
    diff = fabs(diffMatrix.minCoeff()) > diff ? fabs(diffMatrix.minCoeff()) : diff;
    EXPECT_NEAR(0.0, diff, 2.0e-9);
  };
}

} /* namespace Serenity */
