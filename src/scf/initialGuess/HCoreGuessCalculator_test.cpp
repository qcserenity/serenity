/**
 * @file HCoreGuessCalculator_test.cpp
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
#include "scf/initialGuess/HCoreGuessCalculator.h"
#include "data/ElectronicStructure.h"
#include "integrals/OneElectronIntegralController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class HCoreGuessCalculatorTest : public ::testing::Test {
 protected:
  HCoreGuessCalculatorTest() {
  }

  virtual ~HCoreGuessCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(HCoreGuessCalculatorTest, COMinBas) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, true);

  HCoreGuessCalculator<RESTRICTED> guessCalc;

  auto elecStruct = guessCalc.calculateInitialGuess(sys);

  auto& overlaps = elecStruct->getOneElectronIntegralController()->getOverlapIntegrals();

  auto densMat = elecStruct->getDensityMatrix();

  double elecSum = 0;

  for (unsigned int i = 0; i < densMat.rows(); i++) {
    for (unsigned int j = 0; j < densMat.rows(); j++) {
      elecSum += densMat(i, j) * overlaps(i, j);
    }
  }

  EXPECT_NEAR(densMat(0, 0), 2.0978010420755946, 1e-5);
  EXPECT_NEAR(densMat(1, 0), -0.48332073612621285, 1e-5);
  EXPECT_NEAR(densMat(3, 0), -0.38405494649120392, 1e-5);
  EXPECT_NEAR(densMat(5, 0), -0.015089453996145906, 1e-5);
  EXPECT_NEAR(densMat(0, 2), 0.0, 1e-5);
  EXPECT_NEAR(elecSum, 14.0, 1e-5);
}

} // namespace Serenity
