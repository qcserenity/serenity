/**
 * @file AtomicDensityGuessCalculator_test.cpp
 *
 * @date Mar 17, 2017
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
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class AtomicDensityGuessCalculatorTest : public ::testing::Test {
 protected:
  AtomicDensityGuessCalculatorTest() {
  }

  virtual ~AtomicDensityGuessCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(AtomicDensityGuessCalculatorTest, O2MinBasRes_ATOMSCF) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::O2_MINBAS_SING, true);

  AtomicDensityGuessCalculator guessCalc(GUESSMODES::SCF);

  auto elecStruct = guessCalc.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  EXPECT_NEAR(densMat(0, 0), 2.1071094855194836, 1e-6);
  EXPECT_NEAR(densMat(1, 0), -0.46128984211922747, 1e-6);
  EXPECT_NEAR(densMat(3, 3), 0.88024781372049266, 1e-6);
  EXPECT_NEAR(densMat(5, 0), -7.3034439855780155e-05, 1e-6);
  EXPECT_NEAR(densMat(0, 2), 0.0, 1e-6);
}

TEST_F(AtomicDensityGuessCalculatorTest, FRes_INPLACE) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::F_MINUS_6_31Gs, true);

  AtomicDensityGuessCalculator guessCalc(GUESSMODES::SCF_INPLACE);

  auto elecStruct = guessCalc.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  EXPECT_NEAR(densMat(0, 0), 2.0993478760188564, 1e-6);
  EXPECT_NEAR(densMat(1, 0), -0.22378422597398678, 1e-6);
  EXPECT_NEAR(densMat(3, 3), 0.90355455479549251, 1e-6);
  EXPECT_NEAR(densMat(5, 0), 0.0, 1e-6);
  EXPECT_NEAR(densMat(0, 2), -0.27192071672522616, 1e-6);
  // this test actually performs an SCF, and in the process leaves behind a F_FREE system folder inside the
  // TestSystem_F_MINUS_6_31Gs folder
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(sys->getSystemPath() + "F_FREE/", "F_FREE");
}

} // namespace Serenity
