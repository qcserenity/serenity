/**
 * @file   SystemController_test.cpp
 * @author Jan Unsleber
 *
 * @date   Jul 13, 2017
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
#include "system/SystemController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class SystemControllerTest
 * @brief Sets everything up for the tests of SystemController.h/.cpp .
 */
class SystemControllerTest : public ::testing::Test {
 protected:
  SystemControllerTest() {
  }

  virtual ~SystemControllerTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests SystemController.h/.cpp: Tests addition of two systems.
 */
TEST_F(SystemControllerTest, addition_res_res) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  auto result = (*act) + (*env);
  EXPECT_EQ(result->getSpin(), act->getSpin() + env->getSpin());
  EXPECT_EQ(result->getCharge(), act->getCharge() + env->getCharge());
  const auto& cA = act->template getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const auto& cE = env->template getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  Eigen::MatrixXd cR(result->template getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getCoefficients());
  cR.block(0, 0, 4, 1) -= cA.col(0);
  cR.block(4, 1, 4, 1) -= cE.col(0);
  EXPECT_NEAR(0.0, cR.sum(), 1e-12);
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.res.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.res");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.res.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE");
}
/**
 * @test
 * @brief Tests SystemController.h/.cpp: Tests addition of two systems.
 */
TEST_F(SystemControllerTest, addition_res_unres) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  env->setSCFMode(UNRESTRICTED);
  auto result = (*act) + (*env);
  EXPECT_EQ(result->getSpin(), act->getSpin() + env->getSpin());
  EXPECT_EQ(result->getCharge(), act->getCharge() + env->getCharge());
  const auto& cA = act->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const auto& cE = env->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  auto cR(result->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients());
  for_spin(cA, cE, cR) {
    cR_spin.block(0, 0, 4, 1) -= cA_spin.col(0);
    cR_spin.block(4, 1, 4, 1) -= cE_spin.col(0);
    EXPECT_NEAR(0.0, cR_spin.sum(), 1e-12);
  };
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE");
}
/**
 * @test
 * @brief Tests SystemController.h/.cpp: Tests addition of two systems.
 */
TEST_F(SystemControllerTest, addition_unres_res) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  act->setSCFMode(UNRESTRICTED);
  auto result = (*act) + (*env);
  EXPECT_EQ(result->getSpin(), act->getSpin() + env->getSpin());
  EXPECT_EQ(result->getCharge(), act->getCharge() + env->getCharge());
  const auto& cA = act->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const auto& cE = env->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  auto cR(result->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients());
  for_spin(cA, cE, cR) {
    cR_spin.block(0, 0, 4, 1) -= cA_spin.col(0);
    cR_spin.block(4, 1, 4, 1) -= cE_spin.col(0);
    EXPECT_NEAR(0.0, cR_spin.sum(), 1e-12);
  };
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE");
}
/**
 * @test
 * @brief Tests SystemController.h/.cpp: Tests addition of two systems.
 */
TEST_F(SystemControllerTest, addition_unres_unres) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  act->setSCFMode(UNRESTRICTED);
  env->setSCFMode(UNRESTRICTED);
  auto result = (*act) + (*env);
  EXPECT_EQ(result->getSpin(), act->getSpin() + env->getSpin());
  EXPECT_EQ(result->getCharge(), act->getCharge() + env->getCharge());
  const auto& cA = act->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const auto& cE = env->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  auto cR(result->template getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getCoefficients());
  for_spin(cA, cE, cR) {
    cR_spin.block(0, 0, 4, 1) -= cA_spin.col(0);
    cR_spin.block(4, 1, 4, 1) -= cE_spin.col(0);
    EXPECT_NEAR(0.0, cR_spin.sum(), 1e-12);
  };
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/"
              "TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5");
  std::remove("TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE");
}

} /* namespace Serenity */
