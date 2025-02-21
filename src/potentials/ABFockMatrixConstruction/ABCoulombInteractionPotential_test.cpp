/**
 * @file ABCoulombInteractionPotential_test.cpp
 *
 * @date Mai 14, 2018
 * @author: Moritz Bensberg
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
#include "potentials/ABFockMatrixConstruction/ABCoulombInteractionPotential.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "potentials/ERIPotential.h"
#include "potentials/RICoulombPotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class ABCoulombInteractionPotentialTest : public ::testing::Test {
 protected:
  ABCoulombInteractionPotentialTest()
    : systemControllerA(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE)),
      systemControllerB(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE)),
      systemControllerC(SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP)) {
  }

  virtual ~ABCoulombInteractionPotentialTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  /// systems
  std::shared_ptr<SystemController> systemControllerA;
  std::shared_ptr<SystemController> systemControllerB;
  std::shared_ptr<SystemController> systemControllerC;
};

/**
 * @test
 * @brief Tests the restricted construction of the f_AB matrix for the density of C.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_NORI_rFockMatrixABC) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(systemControllerA, basisA, basisB, envDmat);

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_ABC = abCoulIntPot.getMatrix();

  EXPECT_NEAR(f_ABC(0, 0), 1.735096863329e-01, 1e-5);
  EXPECT_NEAR(f_ABC(1, 0), 3.592661749060e-01, 1e-5);
  EXPECT_NEAR(f_ABC(2, 0), 8.443252257335e-01, 1e-5);
  EXPECT_NEAR(f_ABC(0, 1), 4.278766921675e-01, 1e-5);
  EXPECT_NEAR(f_ABC(0, 2), 4.439760900176e-02, 1e-5);
  EXPECT_NEAR(f_ABC(0, 3), 2.591646429591e-01, 1e-5);
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AB matrix for the density of C.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_NORI_uFockMatrixABC) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController());

  ABCoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED> abCoulIntPot(systemControllerA, basisA, basisB, envDmat);

  SPMatrix<Options::SCF_MODES::UNRESTRICTED> f_ABC = abCoulIntPot.getMatrix();

  EXPECT_NEAR(f_ABC.alpha(0, 0), 1.735096863329e-01, 1e-5);
  EXPECT_NEAR(f_ABC.alpha(1, 0), 3.592661749060e-01, 1e-5);
  EXPECT_NEAR(f_ABC.alpha(2, 0), 8.443252257335e-01, 1e-5);
  EXPECT_NEAR(f_ABC.alpha(0, 1), 4.278766921675e-01, 1e-5);
  EXPECT_NEAR(f_ABC.alpha(0, 2), 4.439760900176e-02, 1e-5);
  EXPECT_NEAR(f_ABC.alpha(0, 3), 2.591646429591e-01, 1e-5);

  EXPECT_NEAR(f_ABC.beta(0, 0), 1.735096863329e-01, 1e-5);
  EXPECT_NEAR(f_ABC.beta(1, 0), 3.592661749060e-01, 1e-5);
  EXPECT_NEAR(f_ABC.beta(2, 0), 8.443252257335e-01, 1e-5);
  EXPECT_NEAR(f_ABC.beta(0, 1), 4.278766921675e-01, 1e-5);
  EXPECT_NEAR(f_ABC.beta(0, 2), 4.439760900176e-02, 1e-5);
  EXPECT_NEAR(f_ABC.beta(0, 3), 2.591646429591e-01, 1e-5);
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AA matrix for the density of B.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_NORI_rFockMatrixAAB) {
  auto basisA = systemControllerA->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerB->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(systemControllerA, basisA, basisA, envDmat);

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_AAB = abCoulIntPot.getMatrix();

  EXPECT_NEAR(f_AAB(0, 0), 0.73825203362777958, 1e-5);
  EXPECT_NEAR(f_AAB(1, 0), 0.48242842416584525, 1e-5);
  EXPECT_NEAR(f_AAB(2, 0), 0.44437672997778382, 1e-5);
  EXPECT_NEAR(f_AAB(0, 1), 0.48242842416584525, 1e-5);
  EXPECT_NEAR(f_AAB(0, 2), 0.44437672997778382, 1e-5);
  EXPECT_NEAR(f_AAB(0, 3), 0.40826246586283416, 1e-5);
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AA matrix for the density of A.
 *        Identical to the HFPotential without exchange.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_NORI_rFockMatrixCCC) {
  auto basisC = systemControllerC->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ERIPotential<Options::SCF_MODES::RESTRICTED> hfPotential(
      systemControllerC, envDmat[0], 0.0, systemControllerC->getSettings().basis.integralThreshold, 0.0, 0.0, 0);

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(systemControllerC, basisC, basisC, envDmat);

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_CCCAB = abCoulIntPot.getMatrix();
  SPMatrix<Options::SCF_MODES::RESTRICTED> f_CCC = hfPotential.getMatrix();

  auto diff = (f_CCCAB - f_CCC).array().abs().maxCoeff();
  // only white noise allowed!
  EXPECT_NEAR(diff, 0.0, 1e-12);
}

/**
 * @test
 * @brief Tests the unrestricted construction of the f_AA matrix for the density of A.
 *        Identical to the RICoulombPotential.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_RI_rFockMatrixCCC) {
  auto basisC = systemControllerC->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abRICoulIntPot(
      systemControllerC, basisC, basisC, envDmat, false, Options::DENS_FITS::RI,
      systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB),
      {systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)});

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(systemControllerC, basisC, basisC, envDmat);

  RICoulombPotential<Options::SCF_MODES::RESTRICTED> cccCoulPot(
      systemControllerC, envDmat[0],
      RI_J_IntegralControllerFactory::getInstance().produce(
          systemControllerC->getBasisController(Options::BASIS_PURPOSES::DEFAULT),
          systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
      systemControllerC->getSettings().basis.integralThreshold,
      systemControllerC->getSettings().basis.integralIncrementThresholdStart,
      systemControllerC->getSettings().basis.integralIncrementThresholdEnd,
      systemControllerC->getSettings().basis.incrementalSteps);

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_RICCC = cccCoulPot.getMatrix();

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_RICCC_AB = abRICoulIntPot.getMatrix();
  SPMatrix<Options::SCF_MODES::RESTRICTED> f_CCC_AB = abCoulIntPot.getMatrix();

  Eigen::MatrixXd diff1 = f_RICCC_AB - f_CCC_AB;
  Eigen::MatrixXd diff2 = f_RICCC_AB - f_RICCC;
  EXPECT_NEAR(diff1.array().abs().maxCoeff(), 0.0, 5e-3);
  // only white noise allowed!
  EXPECT_NEAR(diff2.array().abs().maxCoeff(), 0.0, 1e-10);
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AB matrix for the density of C with RI.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_RI_rFockMatrixABC) {
  auto basisA = systemControllerA->getBasisController();
  auto basisB = systemControllerB->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  auto joinedGeometry = std::make_shared<Geometry>();
  *joinedGeometry += *systemControllerA->getGeometry();
  *joinedGeometry += *systemControllerB->getGeometry();
  joinedGeometry->deleteIdenticalAtoms();
  std::shared_ptr<BasisController> actAuxBasis = AtomCenteredBasisControllerFactory::produce(
      joinedGeometry, systemControllerA->getSettings().basis.basisLibPath,
      systemControllerA->getSettings().basis.makeSphericalBasis, false, systemControllerA->getSettings().basis.firstECP,
      systemControllerA->getSettings().basis.auxJLabel);

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(
      systemControllerA, basisA, basisB, envDmat, false, Options::DENS_FITS::RI, actAuxBasis,
      {systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)});

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_ABC = abCoulIntPot.getMatrix();

  EXPECT_NEAR(f_ABC(0, 0), 1.735096863329e-01, 1e-5);
  EXPECT_NEAR(f_ABC(1, 0), 3.592661749060e-01, 1e-5);
  EXPECT_NEAR(f_ABC(2, 0), 8.443252257335e-01, 1e-5);
  EXPECT_NEAR(f_ABC(0, 1), 4.278766921675e-01, 1e-5);
  EXPECT_NEAR(f_ABC(0, 2), 4.439760900176e-02, 1e-5);
  EXPECT_NEAR(f_ABC(0, 3), 2.591646429591e-01, 1e-5);
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AC matrix for the density of C with RI.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_RI_rFockMatrixACC) {
  auto basisA = systemControllerA->getBasisController();
  auto basisC = systemControllerC->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>> envDmat(
      1, systemControllerC->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrixController());

  auto joinedGeometry = std::make_shared<Geometry>();
  *joinedGeometry += *systemControllerA->getGeometry();
  *joinedGeometry += *systemControllerC->getGeometry();
  joinedGeometry->deleteIdenticalAtoms();
  std::shared_ptr<BasisController> actAuxBasis = AtomCenteredBasisControllerFactory::produce(
      joinedGeometry, systemControllerA->getSettings().basis.basisLibPath,
      systemControllerA->getSettings().basis.makeSphericalBasis, false, systemControllerA->getSettings().basis.firstECP,
      systemControllerA->getSettings().basis.auxJLabel);

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot(
      systemControllerA, basisA, basisC, envDmat, false, Options::DENS_FITS::RI, actAuxBasis,
      {systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)});

  ABCoulombInteractionPotential<Options::SCF_MODES::RESTRICTED> abCoulIntPot_NORI(systemControllerA, basisA, basisC, envDmat);

  SPMatrix<Options::SCF_MODES::RESTRICTED> f_ACC = abCoulIntPot.getMatrix();
  SPMatrix<Options::SCF_MODES::RESTRICTED> f_ACC_NORI = abCoulIntPot_NORI.getMatrix();

  EXPECT_NEAR((f_ACC - f_ACC_NORI).sum(), 0.0, 1e-3);

  EXPECT_NEAR(f_ACC(0, 0), 1.486724395567e+00, 1e-5);
  EXPECT_NEAR(f_ACC(1, 0), 0.70729525247340619, 1e-5);
  EXPECT_NEAR(f_ACC(2, 0), 0.56998535470706357, 1e-5);
  EXPECT_NEAR(f_ACC(0, 1), 1.1920521363630558, 1e-5);
  EXPECT_NEAR(f_ACC(0, 2), 0.68014742484020752, 1e-5);
  EXPECT_NEAR(f_ACC(0, 3), 0.0, 1e-5);
}

/**
 * @test
 * @brief Tests the restricted construction of the f_AC matrix for the density of C with RI.
 */
TEST_F(ABCoulombInteractionPotentialTest, H2_dimer_RI_uFockMatrixACC) {
  const auto SPIN = Options::SCF_MODES::UNRESTRICTED;

  auto basisA = systemControllerA->getBasisController();
  auto basisC = systemControllerC->getBasisController();
  std::vector<std::shared_ptr<DensityMatrixController<SPIN>>> envDmat(
      1, systemControllerC->getElectronicStructure<SPIN>()->getDensityMatrixController());

  auto joinedGeometry = std::make_shared<Geometry>();
  *joinedGeometry += *systemControllerA->getGeometry();
  *joinedGeometry += *systemControllerC->getGeometry();
  joinedGeometry->deleteIdenticalAtoms();
  std::shared_ptr<BasisController> actAuxBasis = AtomCenteredBasisControllerFactory::produce(
      joinedGeometry, systemControllerA->getSettings().basis.basisLibPath,
      systemControllerA->getSettings().basis.makeSphericalBasis, false, systemControllerA->getSettings().basis.firstECP,
      systemControllerA->getSettings().basis.auxJLabel);

  ABCoulombInteractionPotential<SPIN> abCoulIntPot(
      systemControllerA, basisA, basisC, envDmat, false, Options::DENS_FITS::RI, actAuxBasis,
      {systemControllerC->getBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)});

  ABCoulombInteractionPotential<SPIN> abCoulIntPot_NORI(systemControllerA, basisA, basisC, envDmat);

  SPMatrix<SPIN> f_ACC = abCoulIntPot.getMatrix();
  SPMatrix<SPIN> f_ACC_NORI = abCoulIntPot_NORI.getMatrix();

  EXPECT_NEAR((f_ACC.alpha - f_ACC_NORI.alpha).sum(), 0.0, 1e-3);
  EXPECT_NEAR((f_ACC.beta - f_ACC_NORI.beta).sum(), 0.0, 1e-3);

  EXPECT_NEAR(f_ACC.alpha(0, 0), 1.486724395567e+00, 1e-5);
  EXPECT_NEAR(f_ACC.alpha(1, 0), 0.70729525247340619, 1e-5);
  EXPECT_NEAR(f_ACC.alpha(2, 0), 0.56998534572991622, 1e-5);
  EXPECT_NEAR(f_ACC.alpha(0, 1), 1.1920521363630558, 1e-5);
  EXPECT_NEAR(f_ACC.alpha(0, 2), 0.68014742484020752, 1e-5);
  EXPECT_NEAR(f_ACC.alpha(0, 3), 0.0, 1e-5);

  EXPECT_NEAR(f_ACC.beta(0, 0), 1.486724395567e+00, 1e-5);
  EXPECT_NEAR(f_ACC.beta(1, 0), 0.70729525247340619, 1e-5);
  EXPECT_NEAR(f_ACC.beta(2, 0), 0.56998534572991622, 1e-5);
  EXPECT_NEAR(f_ACC.beta(0, 1), 1.1920521363630558, 1e-5);
  EXPECT_NEAR(f_ACC.beta(0, 2), 0.68014742484020752, 1e-5);
  EXPECT_NEAR(f_ACC.beta(0, 3), 0.0, 1e-5);
}

} /* namespace Serenity */
