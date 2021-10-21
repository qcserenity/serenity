/**
 * @file CorrespondingOrbitals_test.cpp
 *
 * @date Nov 19, 2020
 * @author Anja Massolle
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
#include "analysis/brokenSymmetry/CorrespondingOrbitals.h"
#include "data/OrbitalController.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class CorrespondingOrbitalsTest : public ::testing::Test {
 protected:
  CorrespondingOrbitalsTest() {
  }

  virtual ~CorrespondingOrbitalsTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(CorrespondingOrbitalsTest, RESTRICTED) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::EthaneA_Def2_SVP_BP86);
  auto coeffBefore = sys->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  LocalizationTask locTask(sys);
  locTask.run();
  auto coeffAfter = sys->getActiveOrbitalController<RESTRICTED>()->getCoefficients();
  auto nOcc = sys->getNOccupiedOrbitals<RESTRICTED>();
  CorrespondingOrbitals<RESTRICTED> corrOrbs(coeffBefore, coeffAfter, nOcc, nOcc);
  EXPECT_NEAR(1.0, corrOrbs.getOverlap(nOcc), 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

TEST_F(CorrespondingOrbitalsTest, UNRESTRICTED) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::MethylRad_def2_SVP_PBE);
  auto coeffBefore = sys->getActiveOrbitalController<UNRESTRICTED>()->getCoefficients();
  LocalizationTask locTask(sys);
  locTask.run();
  auto coeffAfter = sys->getActiveOrbitalController<UNRESTRICTED>()->getCoefficients();
  auto nOcc = sys->getNOccupiedOrbitals<UNRESTRICTED>();
  CorrespondingOrbitals<UNRESTRICTED> corrOrbs(coeffBefore, coeffAfter, nOcc, nOcc);
  auto overlap = corrOrbs.getOverlap(nOcc);
  EXPECT_NEAR(1.0, overlap.alpha, 1e-5);
  EXPECT_NEAR(1.0, overlap.beta, 1e-5);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/