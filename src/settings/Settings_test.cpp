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
#include "dft/functionals/BasicFunctionals.h"
#include "dft/functionals/CompositeFunctionals.h"
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
  EXPECT_EQ(settings.customFunc.basicFunctionals.size(), 0);
}

/**
 * @test
 * @brief Tests print function.
 */
TEST(SettingsTests, Print) {
  Settings settings;
  settings.printSettings();
  EXPECT_EQ(std::remove("default.settings"), 0);
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

/**
 * @test
 * @brief Tests print and load function.
 */
TEST(SettingsTests, PrintAndLoad) {
  Settings settings;
  settings.grid.accuracy = 1;
  settings.printSettings();
  std::ifstream settingsFile;
  settingsFile.open("default.settings", std::ifstream::in);
  Settings loadedSettings(settingsFile);
  EXPECT_EQ(loadedSettings.grid.accuracy, 1);
  EXPECT_EQ(loadedSettings.customFunc.basicFunctionals.size(), 0);
  EXPECT_EQ(std::remove("default.settings"), 0);
}

/**
 * @test
 * @brief Tests print and load function with an active custom functional.
 */
TEST(SettingsTests, PrintAndLoadCustomFunctional) {
  Settings settings;
  settings.grid.accuracy = 1;
  settings.customFunc.mu = 3.14;
  settings.customFunc.basicFunctionals.push_back(BasicFunctionals::BASIC_FUNCTIONALS::X_PBE);
  settings.printSettings();
  std::ifstream settingsFile;
  settingsFile.open("default.settings", std::ifstream::in);
  Settings loadedSettings(settingsFile);
  EXPECT_EQ(loadedSettings.grid.accuracy, 1);
  EXPECT_EQ(loadedSettings.customFunc.basicFunctionals.size(), 1);
  EXPECT_TRUE(loadedSettings.customFunc.basicFunctionals[0] == BasicFunctionals::BASIC_FUNCTIONALS::X_PBE);
  EXPECT_EQ(loadedSettings.customFunc.mu, 3.14);
  EXPECT_EQ(std::remove("default.settings"), 0);
}

} /*namespace Serenity*/
