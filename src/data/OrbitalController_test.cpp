/**
 * @file OrbitalController_test.cpp
 *
 * @date Aug 21, 2018
 * @author Moritz Bensberg
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
#include "data/OrbitalController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "tasks/ScfTask.h"
#include "data/ElectronicStructure.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class OrbitalControllerTest
 * @brief Sets everything up for the tests of OrbitalController.h/.cpp .
 */
class OrbitalControllerTest : public ::testing::Test {
 protected:
  OrbitalControllerTest() {
  }

  virtual ~OrbitalControllerTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};


/**
 * @test
* @brief Tests elimination of linear dependencies from transformation matrix.
 */
TEST_F(OrbitalControllerTest, SCF_WithLinearDependendBasisSet) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto H2Linear = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_LINEARDEPENDENT);

  ScfTask<SPIN> task(H2Linear);
  task.run();

  auto orbs = H2Linear->getElectronicStructure<SPIN>()->getMolecularOrbitals()->getCoefficients();
  auto eps  = H2Linear->getElectronicStructure<SPIN>()->getMolecularOrbitals()->getEigenvalues();

  EXPECT_NEAR(orbs.col(2).sum(), 0.0, 1e-17);
  EXPECT_TRUE(eps(2) > 1.0e10);
}

} /*namespace Serenity*/


