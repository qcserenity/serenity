/**
 * @file InitialGuessCalculator_test.cpp
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
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "data/ElectronicStructure.h"
#include "integrals/OneElectronIntegralController.h"
#include "scf/initialGuess/ExtendedHueckel.h"
#include "scf/initialGuess/HCoreGuessCalculator.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class InitialGuessCalculatorTest : public ::testing::Test {
 protected:
  InitialGuessCalculatorTest() {
  }

  virtual ~InitialGuessCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(InitialGuessCalculatorTest, makeUnresFromResHCore) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, true);

  UnrestrictedFromRestrictedGuess unresGuess(std::make_shared<HCoreGuessCalculator<RESTRICTED>>());

  auto elecStruct = unresGuess.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  double diff = 0;

  for (unsigned int i = 0; i < densMat.alpha.rows(); i++) {
    for (unsigned int j = 0; j < densMat.alpha.rows(); j++) {
      diff += fabs(densMat.alpha(i, j) - densMat.beta(i, j));
    }
  }

  EXPECT_GT(diff, 0.5);
}

TEST_F(InitialGuessCalculatorTest, makeUnresFromResEHT) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS, true);

  UnrestrictedFromRestrictedGuess unresGuess(std::make_shared<ExtendedHueckel>());

  auto elecStruct = unresGuess.calculateInitialGuess(sys);

  auto densMat = elecStruct->getDensityMatrix();

  double diff = 0;

  for (unsigned int i = 0; i < densMat.alpha.rows(); i++) {
    for (unsigned int j = 0; j < densMat.alpha.rows(); j++) {
      diff += fabs(densMat.alpha(i, j) - densMat.beta(i, j));
    }
  }

  EXPECT_GT(diff, 0.5);
}

} // namespace Serenity
