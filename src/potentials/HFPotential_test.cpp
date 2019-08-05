/**
 * @file HFPotential_test.cpp
 *
 * @date Dec 1, 2016
 * @author: Kevin Klahr
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
#include "data/ElectronicStructure.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/HFPotential.h"
#include "potentials/Potential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class HFPotentialTest : public ::testing::Test {
protected:
	HFPotentialTest(): systemController(
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)),
                     systemControllerBP86(
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_BP86)){
	}

	virtual ~HFPotentialTest() = default;

	  static void TearDownTestCase() {
	    SystemController__TEST_SUPPLY::cleanUp();
	  }

	  /// systems
	  std::shared_ptr<SystemController> systemController;
    std::shared_ptr<SystemController> systemControllerBP86;

};

/**
 * @test HFPotentialTest
 * @brief Tests the Fock Matrix of an H2 dimer
 */

TEST_F(HFPotentialTest, H2_FockMatrix) {

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();


  HFPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController,dMat,0.0,0.0);

  FockMatrix<Options::SCF_MODES::RESTRICTED> F = coulPot.getMatrix();

  // TODO create more data for the test
    EXPECT_NEAR(F(0,0), 1.58365844, 1e-5);
    EXPECT_NEAR(F(1,0), 1.038637890, 1e-5);
    EXPECT_NEAR(F(2,0), 0.53449145, 1e-5);
    EXPECT_NEAR(F(0,1), 1.038637890, 1e-5);
    EXPECT_NEAR(F(0,2), 0.53449145, 1e-5);
    EXPECT_NEAR(F(0,3), 0.0, 1e-5);


}

/**
 * @test HFPotentialTest
 * @brief Tests the ERI gradient part of an H2
 */
TEST_F(HFPotentialTest, H2_rgrad) {

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();


  HFPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemController,dMat,1.0,0.0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,2), 0.34677543809435701, 1e-5);
  EXPECT_NEAR(derivative(1,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,2), -0.34677543809435701, 1e-5);

}

/**
 * @test HFPotentialTest
 * @brief Tests the ERI gradient part of an H2
 */
TEST_F(HFPotentialTest, H2_ugrad) {
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemController
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();


  HFPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemController,dMat,1.0,0.0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,2), 0.34677543809435701, 1e-5);
  EXPECT_NEAR(derivative(1,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,2), -0.34677543809435701, 1e-5);

}

/**
 * @test HFPotentialTest
 * @brief Tests the BP86 ERI gradient part of an H2
 */
TEST_F(HFPotentialTest, H2_grad_rBP86) {

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> dMat = systemControllerBP86
      ->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()
      ->getDensityMatrixController();


  HFPotential<Options::SCF_MODES::RESTRICTED> coulPot(systemControllerBP86,dMat,0.0,0.0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,2), 0.71441255550859717, 1e-5);
  EXPECT_NEAR(derivative(1,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,2), -0.71441255550859717, 1e-5);
}

/**
 * @test HFPotentialTest
 * @brief Tests the BP86 ERI gradient part of an H2
 */
TEST_F(HFPotentialTest, H2_grad_uBP86) {

  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMat = systemControllerBP86
      ->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()
      ->getDensityMatrixController();


  HFPotential<Options::SCF_MODES::UNRESTRICTED> coulPot(systemControllerBP86,dMat,0.0,0.0);

  auto derivative = coulPot.getGeomGradients();

  EXPECT_NEAR(derivative(0,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(0,2), 0.71441255550859717, 1e-5);
  EXPECT_NEAR(derivative(1,0), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,1), 0.0, 1e-5);
  EXPECT_NEAR(derivative(1,2), -0.71441255550859717, 1e-5);
}


} /*Namespace*/
