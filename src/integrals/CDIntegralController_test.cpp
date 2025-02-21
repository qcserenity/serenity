/**
 * @file   CDIntegralController_test.cpp
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
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
#include "integrals/CDIntegralController.h"
#include "data/ElectronicStructure.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class CDIntegralController
 * @brief Sets everything up for the tests of CDIntegralController.h/.cpp .
 */
class CDIntegralControllerTest : public ::testing::Test {
 protected:
  CDIntegralControllerTest() {
  }

  virtual ~CDIntegralControllerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Test getter for the CD Threshold
 */
TEST_F(CDIntegralControllerTest, cdThreshold) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  double cdThresh = cdIntController->getDecompositionThreshold();
  EXPECT_NEAR(1e-5, cdThresh, 1e-14);

  SystemController__TEST_SUPPLY::cleanUp();
  cdIntController->cleanup();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Test the diskmode
 */
TEST_F(CDIntegralControllerTest, diskMode) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  EXPECT_FALSE(*cdIntController->getDiskMode());
  cdIntController->setDiskMode();
  EXPECT_TRUE(*cdIntController->getDiskMode());
  cdIntController->unsetDiskMode();
  EXPECT_FALSE(*cdIntController->getDiskMode());

  SystemController__TEST_SUPPLY::cleanUp();

  cdIntController->cleanup();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Generate the ACD basis
 */
TEST_F(CDIntegralControllerTest, generateACD) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  cdIntController->generateACDBasis(systemController->getGeometry());

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());

  cdIntController->cleanup();
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Generate the ACD basis for a larger example
 */
TEST_F(CDIntegralControllerTest, generateACD2) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Ar2_6_31Gs);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  cdIntController->generateACDBasis(systemController->getGeometry());

  std::remove((systemController->getSettings().path + "ACD-6-31GS").c_str());

  cdIntController->cleanup();
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Generate the ACCD basis
 */
TEST_F(CDIntegralControllerTest, generateACCD) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY);
  cdIntController->generateACCDBasis(systemController->getGeometry());

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  std::remove((systemController->getSettings().path + "ACCD-DEF2-TZVP").c_str());

  cdIntController->cleanup();
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: generate the full pseudo ACD vectors
 */
TEST_F(CDIntegralControllerTest, generateACDVectors) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  auto basCont = systemController->getBasisController();
  auto auxBasCont = systemController->getBasisController(Options::BASIS_PURPOSES::ATOMIC_CHOLESKY);
  cdIntController->getACDVectors(basCont, auxBasCont);

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  std::remove((systemController->getSettings().path + "ACCD-DEF2-TZVP").c_str());

  cdIntController->cleanup();
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests CDIntegralController.h/.cpp: Generate pseudo coefficients from density matrix
 */
TEST_F(CDIntegralControllerTest, generatePseudoCoeff) {
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);

  auto cdIntController = systemController->getCDIntegralController();
  assert(cdIntController);

  auto dmatCont = systemController->getElectronicStructure<RESTRICTED>()->getDensityMatrixController();
  auto pseudoCoeff = cdIntController->generatePseudoCoefficients(dmatCont->getDensityMatrix());
  Eigen::MatrixXd dmat = dmatCont->getDensityMatrix();
  auto diff = dmat - pseudoCoeff.second * pseudoCoeff.first.asDiagonal() * pseudoCoeff.second.transpose();

  EXPECT_NEAR(diff.maxCoeff(), 0.0, 1e-14);

  std::remove((systemController->getSettings().path + "ACD-DEF2-TZVP").c_str());
  std::remove((systemController->getSettings().path + "ACCD-DEF2-TZVP").c_str());

  cdIntController->cleanup();
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
