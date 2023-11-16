/**
 * @file   EDATask_test.cpp
 *
 * @date   Aug 22, 2017
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
#include "tasks/EDATask.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class EDATaskTest
 * @brief Sets the EDATask test up.
 */
class EDATaskTest : public ::testing::Test {
 protected:
  EDATaskTest()
    : _systemA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_6_31Gs)),
      _systemB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_6_31Gs)) {
    _path = "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs/";
  }
  virtual ~EDATaskTest() {
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.energies.res").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.orbs.res.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.dmat.res.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.energies.unres").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.orbs.unres.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.dmat.unres.h5").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.settings").c_str());
    std::remove((_path + "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs.xyz").c_str());
    std::remove((_path).c_str());
  }
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /* The systems. */
  std::shared_ptr<SystemController> _systemA;
  std::shared_ptr<SystemController> _systemB;
  std::string _path;
};

/**
 * @test
 * @brief Tests the EDA.
 *
 * The calculations have been verified by the equivalent GAMESS calculations and are in range of the
 * results obtained by Morokuma and Kitaura.
 * (Kitaura-Morokuma : K. Kitaura and K. Morokuma, Int. J. Quantum Chem. 10, 325 (1976))
 *
 */
TEST_F(EDATaskTest, EDA_rWaterDimer) {
  const auto SCFMode = Options::SCF_MODES::RESTRICTED;
  EDATask<SCFMode> edaTask({_systemA, _systemB});
  edaTask.run();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(_path,
                                                        "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs");
}
/**
 * @test
 * @brief Tests the EDA.
 *
 * The calculations have been verified by the equivalent GAMESS calculations and are in range of the
 * results obtained by Morokuma and Kitaura.
 * (Kitaura-Morokuma : K. Kitaura and K. Morokuma, Int. J. Quantum Chem. 10, 325 (1976))
 *
 */
TEST_F(EDATaskTest, EDA_uWaterDimer) {
  const auto SCFMode = Options::SCF_MODES::UNRESTRICTED;

  EDATask<SCFMode> edaTask({_systemA, _systemB});
  edaTask.run();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(_path,
                                                        "TestSystem_WaterMonOne_6_31Gs+TestSystem_WaterMonTwo_6_31Gs");
}
} /* namespace Serenity */
