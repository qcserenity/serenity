/**
 * @file ExternalChargeController_test.cpp
 *
 * @date Apr. 29, 2024
 * @author Moritz Bensberg
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
#include "data/ExternalChargeController.h"
#include "geometry/Point.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ExternalChargeControllerTest : public ::testing::Test {
 protected:
  ExternalChargeControllerTest() {
  }

  virtual ~ExternalChargeControllerTest() = default;
};

/**
 * @brief Test the external charge controller.
 */
TEST_F(ExternalChargeControllerTest, readChargesFromFile) {
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/TestSystem_WCCR1010_P1_Def2_SVP_HF/external_charges.dat";
  }
  else {
    throw SerenityError("ERROR: Environment variable SERENITY_RESOURCES not set.");
  }

  Settings settings;
  settings.extCharges.externalChargesFile = pathToTestsResources;
  ExternalChargeController chargeController(settings);
  const std::vector<std::pair<double, Point>> charges = chargeController.getExternalCharges();
  EXPECT_EQ(charges.size(), 4803);
  EXPECT_NEAR(charges[0].first, -0.573704, 5e-6);
  EXPECT_NEAR(charges[1].first, 0.290609, 5e-6);
}

} /* namespace Serenity */
