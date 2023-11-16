/**
 * @file ABERIPotential_test.cpp
 *
 * @date Jun 23, 2018
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
#include "potentials/ABFockMatrixConstruction/ABERIPotential.h"
#include "data/ElectronicStructure.h"
#include "potentials/ERIPotential.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABERIPotentialTest : public ::testing::Test {
 protected:
  ABERIPotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~ABERIPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
};

/**
 * @test
 * @brief Tests the restricted construction of the f_AA matrix.
 */
TEST_F(ABERIPotentialTest, H2_dimer_rFockMatrixAA) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto auxBasisA = systemControllerA->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> densityMatrixController = {
      systemControllerA->getElectronicStructure<SPIN>()->getDensityMatrixController()};

  ABERIPotential<SPIN> abHFPot(systemControllerA, basisA, basisA, densityMatrixController, 1.0, 0.0, 0.0, false,
                               Options::DENS_FITS::RI, auxBasisA);

  ERIPotential<SPIN> hfPotential(systemControllerA,
                                 systemControllerA->getElectronicStructure<SPIN>()->getDensityMatrixController(), 1.0,
                                 1E-10, 0.0, 0.0, 0, false);

  SPMatrix<SPIN> f_AA1 = abHFPot.getMatrix();
  SPMatrix<SPIN> f_AA2 = hfPotential.getMatrix();

  SPMatrix<SPIN> test = f_AA1.transpose();

  EXPECT_NEAR((f_AA1 - test).sum(), 0.0, 1e-10);
  EXPECT_NEAR(f_AA1(0, 0), f_AA2(0, 0), 1e-5);
  EXPECT_NEAR(f_AA1(1, 0), f_AA2(1, 0), 1e-5);
  EXPECT_NEAR(f_AA1(2, 0), f_AA2(2, 0), 1e-5);
  EXPECT_NEAR(f_AA1(0, 1), f_AA2(0, 1), 1e-5);
  EXPECT_NEAR(f_AA1(0, 2), f_AA2(0, 2), 1e-5);
  EXPECT_NEAR(f_AA1(0, 3), f_AA2(0, 3), 1e-5);
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AB matrix.
 */
TEST_F(ABERIPotentialTest, H2_dimer_rFockMatrixAB) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> densityMatrixController = {
      systemControllerA->getElectronicStructure<SPIN>()->getDensityMatrixController()};

  ABERIPotential<SPIN> abHFPot(systemControllerA, basisA, basisB, densityMatrixController, 1.0);

  SPMatrix<SPIN> f_AA1 = abHFPot.getMatrix();

  EXPECT_NEAR(f_AA1(0, 0), 0.96824147145617268, 1e-5);
  EXPECT_NEAR(f_AA1(1, 0), 0.2964919220058731, 1e-5);
  EXPECT_NEAR(f_AA1(2, 0), 0.17671512858234784, 1e-5);
  EXPECT_NEAR(f_AA1(0, 1), 0.64405197133753755, 1e-5);
  EXPECT_NEAR(f_AA1(0, 2), 0.28593915743733905, 1e-5);
  EXPECT_NEAR(f_AA1(0, 3), 0.0, 1e-5);
}

} /* namespace Serenity */
