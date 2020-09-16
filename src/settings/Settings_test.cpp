/**
 * @file Settings_test.cpp
 *
 * @date Aug 22, 2017
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
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <istream>

namespace Serenity {
/**
 * @test
 * @brief Tests reading settings from file (stream);
 */
TEST(SettingsTests, ReadFromFile) {
  std::string pathToTestsResources;
  if (const char* env_p = std::getenv("SERENITY_RESOURCES")) {
    pathToTestsResources = (std::string)env_p + "testresources/";
  }
  else {
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set." << std::endl;
    assert(false);
  }
  std::ifstream newSettingsFile;
  newSettingsFile.open(pathToTestsResources + "TestSystem_H2_6_31Gs_BP86/TestSystem_H2_6_31Gs_BP86.settings",
                       std::ifstream::in);
  Settings settings(newSettingsFile);
  EXPECT_EQ("TestSystem_H2_6_31Gs_BP86", settings.name);
  EXPECT_TRUE(settings.dft.functional == CompositeFunctionals::XCFUNCTIONALS::BP86);
}

/**
 * @test
 * @brief Tests print function.
 */
TEST(SettingsTests, Print) {
  Settings settings;
  settings.printSettings();
  std::remove("default.settings");
}
/**
 * @test
 * @brief Tests set function using all strings.
 */
TEST(SettingsTests, Set) {
  Settings settings;
  settings.set("DFT", "functional", "pbe0");
  EXPECT_TRUE(settings.dft.functional == CompositeFunctionals::XCFUNCTIONALS::PBE0);
}

} /*namespace Serenity*/
