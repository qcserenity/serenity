/**
 * @file SuperpositionOfAtomicPotentials_test.cpp
 *
 * @date September 8, 2019
 * @author Jan Unsleber
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
#include "scf/initialGuess/SuperpositionOfAtomicPotentials.h"
#include "data/ElectronicStructure.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class SuperpositionOfAtomicPotentialsTest : public ::testing::Test {
 protected:
  SuperpositionOfAtomicPotentialsTest() {
  }

  virtual ~SuperpositionOfAtomicPotentialsTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SuperpositionOfAtomicPotentialsTest, O2MinBasRes) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING, true);

  SuperpositionOfAtomicPotentials guessCalc;

  auto elecStruct = guessCalc.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  EXPECT_NEAR(densMat(0, 0), +2.1110200656169869, 1e-5);
  EXPECT_NEAR(densMat(1, 0), -0.4847541905973736, 1e-5);
  EXPECT_NEAR(densMat(3, 3), +0.8363779864504210, 1e-5);
  EXPECT_NEAR(densMat(5, 0), -0.0039836145553974, 1e-5);
  EXPECT_NEAR(densMat(0, 2), 0.0, 1e-5);
}

TEST_F(SuperpositionOfAtomicPotentialsTest, FRes) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs, true);

  SuperpositionOfAtomicPotentials guessCalc;

  auto elecStruct = guessCalc.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  EXPECT_NEAR(densMat(0, 0), +2.0967977919552143, 1e-5);
  EXPECT_NEAR(densMat(1, 0), -0.2073679433347872, 1e-5);
  EXPECT_NEAR(densMat(3, 3), +0.9459396557309102, 1e-5);
  EXPECT_NEAR(densMat(5, 0), 0.0, 1e-5);
  EXPECT_NEAR(densMat(0, 2), -0.2858770514488218, 1e-5);
}

} // namespace Serenity
