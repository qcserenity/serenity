/**
 * @file ELFCalculator_test.cpp
 *
 * @date Apr 19, 2016, rework Nov 23, 2018
 * @author Kevin Klahr, Melanie Börner
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
#include "analysis/localizationFunctions/ELFCalculator.h"
#include "data/grid/GridData.h"
#include "grid/GridController.h"
#include "system/SystemController.h"
#include "tasks/ScfTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
template<Options::SCF_MODES SPIN>
class ELFCalculator;

class ELFCalculatorTest : public ::testing::Test {
 protected:
  ELFCalculatorTest()
    : systemController(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS)) {
  }
  virtual ~ELFCalculatorTest() = default;
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
  /// system
  std::shared_ptr<SystemController> systemController;
};

/**
 * @test
 * @brief Tests the calculation of the ELF of CO in a minimal basis
 */
TEST_F(ELFCalculatorTest, TotalELF) {
  ELFCalculator<Options::SCF_MODES::RESTRICTED> calcres(systemController);
  ELFCalculator<Options::SCF_MODES::UNRESTRICTED> calcunres(systemController);

  auto grid = systemController->getGridController(Options::GRID_PURPOSES::SMALL);

  auto elfRes = calcres.calculateTotalELFOnGrid(grid);
  auto elfUnres = calcunres.calculateTotalELFOnGrid(grid);

  for (unsigned int i = 0; i < elfRes.size(); i++) {
    EXPECT_NEAR(elfRes[i], elfUnres[i], 1e-10);
  }
}
/**
 * @test
 * @brief Tests the calculation of the ELF of CO in a minimal basis
 */
TEST_F(ELFCalculatorTest, TotalELFTS) {
  ELFCalculator<Options::SCF_MODES::RESTRICTED> calcres(systemController);
  ELFCalculator<Options::SCF_MODES::UNRESTRICTED> calcunres(systemController);

  auto grid = systemController->getGridController(Options::GRID_PURPOSES::SMALL);

  auto elfRes = calcres.calculateTotalELFTSOnGrid(grid);
  auto elfUnres = calcunres.calculateTotalELFTSOnGrid(grid);

  for (unsigned int i = 0; i < elfRes.size(); i++) {
    // less accurate due to 2nd derivatives of density on grid
    EXPECT_NEAR(elfRes[i], elfUnres[i], 1e-6);
  }
}

/*
 * The reference values for the following tests have been calculated with:
 *
 * import numpy as np
 *
 * {Declare Tau, Rho, gradRho2 and grad2Rho}
 *
 * #ELF
 * D_p=Tau-((1.0*gradRho2)/(8.0*Rho))
 * D_p0=(3.0/10.0)*((3.0*np.pi*np.pi)**(2.0/3.0))*(Rho**(5.0/3.0))
 * ELF=1.0/(1.0+((D_p/D_p0)**2.0))
 * print(ELF)
 *
 * #ELF TS
 * Tau_k=(3.0/10.0)*((3.0*np.pi*np.pi)**(2.0/3.0))*(Rho**(5.0/3.0))+((1.0*gradRho2)/(72.0*Rho))+(grad2Rho/6.0)
 * D_TS=Tau_k-((1.0*gradRho2)/(8.0*Rho))
 * ELF_TS=1.0/(1.0+((D_TS/D_p0)**2.0))
 * print(ELF_TS)
 *
 */

/**
 * @test
 * @brief Tests specific gridpoints of the ELF of CO in a minimal basis
 * Note that the mock data to create these reference values was obtained in parallel on 16 cores,
 * which, unfortunately, makes a difference due to grid parallelization. This means that the test also has to
 * run on 16 cores or might otherwise fail in very random fashions. As i am writing this, the GitLab
 * server runs all tests on 16 cores, but since this might change in the future, heed my warning.
 */
// TEST_F(ELFCalculatorTest, TotalELFReference) {
//  ELFCalculator<Options::SCF_MODES::RESTRICTED> calcres(systemController);
//
//  auto grid = systemController->getGridController(Options::GRID_PURPOSES::SMALL);
//
//  auto elfRes= calcres.calculateTotalELFOnGrid(grid);
//
//  EXPECT_NEAR(elfRes[100],0.0158037405615,1e-5);
//  EXPECT_NEAR(elfRes[200],0.869145707939,1e-5);
//  EXPECT_NEAR(elfRes[300],0.651072177177,1e-5);
//  EXPECT_NEAR(elfRes[400],0.00616933378342,1e-5);
//
//}

/**
 * @test
 * @brief Tests specific gridpoints  of the ELF TS of CO in a minimal basis
 * Note that the mock data to create these reference values was obtained in parallel on 16 cores,
 * which, unfortunately, makes a difference due to grid parallelization. This means that the test also has to
 * run on 16 cores or might otherwise fail in very random fashions. As i am writing this, the GitLab
 * server runs all tests on 16 cores, but since this might change in the future, heed my warning.
 */
// TEST_F(ELFCalculatorTest, TotalELFTSReference) {
//  ELFCalculator<Options::SCF_MODES::RESTRICTED> calcres(systemController);
//
//  auto grid = systemController->getGridController(Options::GRID_PURPOSES::SMALL);
//
//  auto elfRes= calcres.calculateTotalELFTSOnGrid(grid);
//
//  EXPECT_NEAR(elfRes[100],0.00099875752324,1e-5);
//  EXPECT_NEAR(elfRes[200],0.790986787932,1e-5);
//  EXPECT_NEAR(elfRes[300],0.510729219829,1e-5);
//  EXPECT_NEAR(elfRes[400],0.00015923561279,1e-5);
//
//}

} // namespace Serenity
