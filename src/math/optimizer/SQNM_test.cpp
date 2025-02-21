/**
 * @file   SQNM_test.cpp
 *
 * @date   Oct 09, 2024
 * @author Thorben Wiegmann
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
#include "math/optimizer/SQNM.h"
#include "data/ElectronicStructure.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/GeometryOptimizationTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

// All results were obtained using pysisyphus (Version 0.7.6.post2) for the optimization and ORCA (Version 6.0.0) for
// the energy and gradient calculations.

namespace Serenity {

class SQNMTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SQNMTest, SQNMConvergence) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Water_Dimer_def2_SVP_HF, true);
  std::vector<std::shared_ptr<SystemController>> activeSystems = {system};
  std::vector<std::shared_ptr<SystemController>> passiveSystems = {};
  GeometryOptimizationTask<Options::SCF_MODES::RESTRICTED> task(activeSystems, passiveSystems);
  task.settings.optAlgorithm = Options::OPTIMIZATION_ALGORITHMS::SQNM;
  task.run();
  double energy = system->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy();
  // Keep in mind that pysisyphus uses different convergence criteria.
  EXPECT_NEAR(energy, -151.93197739, 1e-5);
  EXPECT_EQ(0, std::remove((system->getSystemPath() + system->getSystemName() + ".trj").c_str()));
}
} // namespace Serenity
