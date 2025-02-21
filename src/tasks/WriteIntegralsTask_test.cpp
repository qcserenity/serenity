/**
 * @file WriteIntegralsTask_test.cpp
 *
 * @date   Nov 10, 2023
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
#include "tasks/WriteIntegralsTask.h"                 //To be tested.
#include "system/SystemController.h"                  //Get system name.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.
#include <string>        //file name

namespace Serenity {
class WriteIntegralsTaskTest : public ::testing::Test {
 protected:
  WriteIntegralsTaskTest() {
  }
  virtual ~WriteIntegralsTaskTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(WriteIntegralsTaskTest, writeHCoreIntegralsHDF5) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  WriteIntegralsTask writeIntegralsTask(sys);
  writeIntegralsTask.settings.hCoreIntegrals = true;
  writeIntegralsTask.settings.fileFormat = Options::INTEGRAL_FILE_TYPES::HDF5;
  writeIntegralsTask.run();

  std::string fileName = sys->getSystemName() + ".hcore.h5";
  std::ifstream input(fileName.c_str());
  ASSERT_TRUE(input.good());
  std::remove(fileName.c_str());
}

TEST_F(WriteIntegralsTaskTest, writeHCoreIntegralsASCII) {
  auto sys = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WCCR1010_def2_SVP_HF);
  WriteIntegralsTask writeIntegralsTask(sys);
  writeIntegralsTask.settings.hCoreIntegrals = true;
  writeIntegralsTask.settings.fileFormat = Options::INTEGRAL_FILE_TYPES::ASCII;
  writeIntegralsTask.run();

  std::string fileName = sys->getSystemName() + ".hcore.dat";
  std::ifstream input(fileName.c_str());
  ASSERT_TRUE(input.good());
  std::remove(fileName.c_str());
}

} /*namespace Serenity*/
