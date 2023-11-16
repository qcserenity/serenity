/**
 * @file FosterBoysLocalization_test.cpp
 *
 * @date Dec 2, 2015
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
#include "analysis/orbitalLocalization/FosterBoysLocalization.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/Matrix.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <iostream>

namespace Serenity {

class FosterBoysLocalizationTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test FB Localization Restricted CO minimal basis
 * TODO check whether the density matrix stays the same
 */
TEST_F(FosterBoysLocalizationTest, testLocalizationRestricted) {
  // Create the test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  // Perform SCF
  system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  // Create a copy of the SCF orbitals
  auto orbitals = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(
      *system->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>());
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients = orbitals->getCoefficients();

  const auto& oneIntController = system->getOneElectronIntegralController();
  const auto& overlaps = oneIntController->getOverlapIntegrals();
  /*
   * Calculate a simple Pipek-Mezey-like localization measure to check whether the algorithm does
   * localize something.
   */
  double convMeasure = 0.0;
  for (unsigned int i = 0; i < 7; ++i) {
    double sumCoeffsThisOrbFirstAtom = 0.0;
    for (unsigned int mu = 0; mu < 5; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        sumCoeffsThisOrbFirstAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasure += sumCoeffsThisOrbFirstAtom * sumCoeffsThisOrbFirstAtom;
    double sumCoeffsThisOrbSecondAtom = 0.0;
    for (unsigned int mu = 5; mu < 10; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        sumCoeffsThisOrbSecondAtom += coefficients(mu, i) * coefficients(nu, i) * overlaps(mu, nu);
      }
    }
    convMeasure += sumCoeffsThisOrbSecondAtom * sumCoeffsThisOrbSecondAtom;
  }
  FosterBoysLocalization<Options::SCF_MODES::RESTRICTED> localizer(system);
  /*
   * Perform localization
   */
  auto nOcc = system->getNOccupiedOrbitals<RESTRICTED>();
  std::vector<unsigned int> range;
  for (unsigned int i = 0; i < nOcc; ++i)
    range.push_back(i);
  localizer.localizeOrbitals(*orbitals, 500, range);
  coefficients = orbitals->getCoefficients();
  /*
   * check if MOs still orthogonal
   */
  Eigen::MatrixXd overlapMatrix = coefficients.transpose() * overlaps * coefficients;
  Eigen::MatrixXd identity(orbitals->getNOrbitals(), orbitals->getNOrbitals());
  identity.setIdentity();
  Eigen::MatrixXd diffMatrix = overlapMatrix - identity;
  double diff = fabs(diffMatrix.maxCoeff());
  diff = fabs(diffMatrix.minCoeff()) > diff ? fabs(diffMatrix.minCoeff()) : diff;
  EXPECT_NEAR(0.0, diff, 2.0e-9);
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
  EXPECT_GT(convMeasureLoc, convMeasure + 0.1);
}

// Test FB Localization Unrestricted
TEST_F(FosterBoysLocalizationTest, testLocalizationUnrestricted) {
  // Create the test system
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  // Perform SCF
  system->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  // Create a copy of the SCF orbitals
  auto orbitals = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(
      *system->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>());
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficients = orbitals->getCoefficients();

  const auto& oneIntController = system->getOneElectronIntegralController();
  const auto& overlaps = oneIntController->getOverlapIntegrals();
  /*
   * Calculate a simple Pipek-Mezey-like localization measure to check whether the algorithm does
   * localize something.
   */
  double convMeasure = 0.0;
  for (unsigned int i = 0; i < 7; ++i) {
    double sumCoeffsThisOrbFirstAtom = 0.0;
    for (unsigned int mu = 0; mu < 5; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        for_spin(coefficients) {
          sumCoeffsThisOrbFirstAtom += coefficients_spin(mu, i) * coefficients_spin(nu, i) * overlaps(mu, nu);
        };
      }
    }
    convMeasure += sumCoeffsThisOrbFirstAtom * sumCoeffsThisOrbFirstAtom;
    double sumCoeffsThisOrbSecondAtom = 0.0;
    for (unsigned int mu = 5; mu < 10; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        for_spin(coefficients) {
          sumCoeffsThisOrbSecondAtom += coefficients_spin(mu, i) * coefficients_spin(nu, i) * overlaps(mu, nu);
        };
      }
    }
    convMeasure += sumCoeffsThisOrbSecondAtom * sumCoeffsThisOrbSecondAtom;
  }

  FosterBoysLocalization<Options::SCF_MODES::UNRESTRICTED> localizer(system);
  /*
   * Perform localization
   */
  auto nOcc = system->getNOccupiedOrbitals<UNRESTRICTED>();
  SpinPolarizedData<UNRESTRICTED, std::vector<unsigned int>> range;
  for_spin(range, nOcc) {
    for (unsigned int i = 0; i < nOcc_spin; ++i)
      range_spin.push_back(i);
  };
  localizer.localizeOrbitals(*orbitals, 500, range);
  coefficients = orbitals->getCoefficients();
  /*
   * check if MOs still orthogonal
   */
  for_spin(coefficients) {
    Eigen::MatrixXd overlapMatrix = coefficients_spin.transpose() * overlaps * coefficients_spin;
    Eigen::MatrixXd identity(orbitals->getNOrbitals(), orbitals->getNOrbitals());
    identity.setIdentity();
    Eigen::MatrixXd diffMatrix = overlapMatrix - identity;
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
        for_spin(coefficients) {
          sumCoeffsThisOrbFirstAtom += coefficients_spin(mu, i) * coefficients_spin(nu, i) * overlaps(mu, nu);
        };
      }
    }
    convMeasureLoc += sumCoeffsThisOrbFirstAtom * sumCoeffsThisOrbFirstAtom;
    double sumCoeffsThisOrbSecondAtom = 0.0;
    for (unsigned int mu = 5; mu < 10; ++mu) {
      for (unsigned int nu = 0; nu < 10; ++nu) {
        for_spin(coefficients) {
          sumCoeffsThisOrbSecondAtom += coefficients_spin(mu, i) * coefficients_spin(nu, i) * overlaps(mu, nu);
        };
      }
    }
    convMeasureLoc += sumCoeffsThisOrbSecondAtom * sumCoeffsThisOrbSecondAtom;
  }
  EXPECT_GT(convMeasureLoc, convMeasure + 0.1);
}

} /* namespace Serenity */
