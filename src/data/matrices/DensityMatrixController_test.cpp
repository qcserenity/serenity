/**
 * @file DensityMatrixController_test.cpp
 *
 * @date May 31, 2016
 * @author: Jan Unsleber
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
#include "data/matrices/DensityMatrixController.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class DensityMatrixControllerNotificationTest : public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  DensityMatrixControllerNotificationTest() {
    notified = false;
  }
  ~DensityMatrixControllerNotificationTest() = default;
  void notify() {
    notified = true;
  }
  bool notified;

 protected:
};

class DensityMatrixControllerTest : public ::testing::Test {
 protected:
  DensityMatrixControllerTest()
    : tBasisController(BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::MINIMAL)),
      tRestrictedDensityMatrix(tBasisController),
      tUnrestrictedDensityMatrix(tBasisController) {
    tRestrictedDensityMatrix(0, 0) = 0.08;
    tRestrictedDensityMatrix(0, 1) = 0.32;
    tRestrictedDensityMatrix(1, 0) = 0.32;
    tRestrictedDensityMatrix(1, 1) = 1.28;
    tUnrestrictedDensityMatrix.alpha(0, 0) = 0.04;
    tUnrestrictedDensityMatrix.alpha(0, 1) = 0.16;
    tUnrestrictedDensityMatrix.alpha(1, 0) = 0.16;
    tUnrestrictedDensityMatrix.alpha(1, 1) = 0.64;
    tUnrestrictedDensityMatrix.beta(0, 0) = 0.04;
    tUnrestrictedDensityMatrix.beta(0, 1) = 0.16;
    tUnrestrictedDensityMatrix.beta(1, 0) = 0.16;
    tUnrestrictedDensityMatrix.beta(1, 1) = 0.64;
    restrictedOccupations = Eigen::VectorXd(2);
    restrictedOccupations[0] = 2.0;
    restrictedOccupations[1] = 0.0;
    unrestrictedOccupations = makeUnrestrictedFromPieces<Eigen::VectorXd>(Eigen::VectorXd(2), Eigen::VectorXd(2));
    unrestrictedOccupations.alpha[0] = 1.0;
    unrestrictedOccupations.alpha[1] = 0.0;
    unrestrictedOccupations.beta[0] = 1.0;
    unrestrictedOccupations.beta[1] = 0.0;
  }
  virtual ~DensityMatrixControllerTest() = default;
  // The basis in which it is defined
  std::shared_ptr<BasisController> tBasisController;
  // The object under test
  DensityMatrix<Options::SCF_MODES::RESTRICTED> tRestrictedDensityMatrix;
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> tUnrestrictedDensityMatrix;
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd> restrictedOccupations;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> unrestrictedOccupations;
};
//=====================
//     RESTRICTED
//=====================
/**
 * @test
 * @brief Tests DensityMatrix.h: Construction from DensityMatrix (restricted)
 */
TEST_F(DensityMatrixControllerTest, Constructor_DensityMatrix_R) {
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(tRestrictedDensityMatrix);
  auto dmat = controller.getDensityMatrix();
  EXPECT_EQ(0.08, dmat(0, 0));
  EXPECT_EQ(0.32, dmat(1, 0));
  EXPECT_EQ(0.32, dmat(0, 1));
  EXPECT_EQ(1.28, dmat(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Construction from Orbitals (restricted)
 */
TEST_F(DensityMatrixControllerTest, Constructor_Orbitals_R) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::RESTRICTED>>();
  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients = orbset->getCoefficients();
  coefficients(0, 0) = 0.2;
  coefficients(0, 1) = 0.9;
  coefficients(1, 0) = 0.8;
  coefficients(1, 1) = 0.1;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  EXPECT_FALSE(nTest->notified);
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(orbset, restrictedOccupations);
  controller.addSensitiveObject(nTest);
  auto dmat = controller.getDensityMatrix();
  EXPECT_DOUBLE_EQ(0.08, dmat(0, 0));
  EXPECT_DOUBLE_EQ(0.32, dmat(1, 0));
  EXPECT_DOUBLE_EQ(0.32, dmat(0, 1));
  EXPECT_DOUBLE_EQ(1.28, dmat(1, 1));
  controller.updateDensityMatrix();
  EXPECT_FALSE(nTest->notified);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Attach an OrbitalController (restricted)
 */
TEST_F(DensityMatrixControllerTest, AttachOrbitalController_OccVec_R) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::RESTRICTED>>();
  EXPECT_FALSE(nTest->notified);
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(tRestrictedDensityMatrix);
  controller.addSensitiveObject(nTest);

  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients = orbset->getCoefficients();
  coefficients(0, 0) = 0.1;
  coefficients(1, 0) = 0.9;
  coefficients(0, 1) = 0.8;
  coefficients(1, 1) = 0.2;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());

  controller.attachOrbitals(orbset, restrictedOccupations, true);

  auto dmat = controller.getDensityMatrix();
  EXPECT_FALSE(nTest->notified);
  EXPECT_DOUBLE_EQ(0.02, dmat(0, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(1, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(0, 1));
  EXPECT_DOUBLE_EQ(1.62, dmat(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Attach an OrbitalController (restricted)
 */
TEST_F(DensityMatrixControllerTest, AttachOrbitalController_OccOrbNumber_R) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::RESTRICTED>>();
  EXPECT_FALSE(nTest->notified);
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(tRestrictedDensityMatrix);
  controller.addSensitiveObject(nTest);

  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients = orbset->getCoefficients();
  coefficients(0, 0) = 0.1;
  coefficients(1, 0) = 0.9;
  coefficients(0, 1) = 0.8;
  coefficients(1, 1) = 0.2;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());

  controller.attachOrbitals(orbset, SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>(1), true);

  auto dmat = controller.getDensityMatrix();
  EXPECT_FALSE(nTest->notified);
  EXPECT_DOUBLE_EQ(0.02, dmat(0, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(1, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(0, 1));
  EXPECT_DOUBLE_EQ(1.62, dmat(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Update via setter (restricted)
 */
TEST_F(DensityMatrixControllerTest, MatrixSetter_R) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::RESTRICTED>>();
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(tRestrictedDensityMatrix);
  controller.addSensitiveObject(nTest);

  DensityMatrix<Options::SCF_MODES::RESTRICTED> newMat(tBasisController);
  newMat(0, 0) = 1.08;
  newMat(0, 1) = 2.32;
  newMat(1, 0) = 2.32;
  newMat(1, 1) = 4.28;

  EXPECT_FALSE(nTest->notified);
  controller.setDensityMatrix(newMat);
  auto dmat = controller.getDensityMatrix();
  EXPECT_TRUE(nTest->notified);
  EXPECT_EQ(1.08, dmat(0, 0));
  EXPECT_EQ(2.32, dmat(1, 0));
  EXPECT_EQ(2.32, dmat(0, 1));
  EXPECT_EQ(4.28, dmat(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Update via notification (restricted)
 */
TEST_F(DensityMatrixControllerTest, NotifiedUpdate_R) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::RESTRICTED>>();
  EXPECT_FALSE(nTest->notified);
  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::RESTRICTED>>(tBasisController, 0);
  // set initial coefficients
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients = orbset->getCoefficients();
  coefficients(0, 0) = 0.2;
  coefficients(0, 1) = 0.9;
  coefficients(1, 0) = 0.8;
  coefficients(1, 1) = 0.1;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  DensityMatrixController<Options::SCF_MODES::RESTRICTED> controller(orbset, restrictedOccupations);
  controller.addSensitiveObject(nTest);

  // check initial
  auto dmat = controller.getDensityMatrix();
  EXPECT_DOUBLE_EQ(0.08, dmat(0, 0));
  EXPECT_DOUBLE_EQ(0.32, dmat(1, 0));
  EXPECT_DOUBLE_EQ(0.32, dmat(0, 1));
  EXPECT_DOUBLE_EQ(1.28, dmat(1, 1));

  // edit obset
  coefficients(0, 0) = 0.1;
  coefficients(1, 0) = 0.9;
  coefficients(0, 1) = 0.8;
  coefficients(1, 1) = 0.2;
  EXPECT_FALSE(nTest->notified);

  // trigger orbset notification by setting new coefficients
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  dmat = controller.getDensityMatrix();
  EXPECT_TRUE(nTest->notified);
  EXPECT_DOUBLE_EQ(0.02, dmat(0, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(1, 0));
  EXPECT_DOUBLE_EQ(0.18, dmat(0, 1));
  EXPECT_DOUBLE_EQ(1.62, dmat(1, 1));
}
//=====================
//    UNRESTRICTED
//=====================
/**
 * @test
 * @brief Tests DensityMatrix.h: Construction from DensityMatrix (unrestricted)
 */
TEST_F(DensityMatrixControllerTest, Constructor_DensityMatrix_U) {
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> controller(tUnrestrictedDensityMatrix);
  auto dmat = controller.getDensityMatrix();
  EXPECT_EQ(0.04, dmat.alpha(0, 0));
  EXPECT_EQ(0.16, dmat.alpha(1, 0));
  EXPECT_EQ(0.16, dmat.alpha(0, 1));
  EXPECT_EQ(0.64, dmat.alpha(1, 1));
  EXPECT_EQ(0.04, dmat.beta(0, 0));
  EXPECT_EQ(0.16, dmat.beta(1, 0));
  EXPECT_EQ(0.16, dmat.beta(0, 1));
  EXPECT_EQ(0.64, dmat.beta(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Construction from Orbitals (unrestricted)
 */
TEST_F(DensityMatrixControllerTest, Constructor_Orbitals_U) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::UNRESTRICTED>>();
  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficients = orbset->getCoefficients();
  coefficients.alpha(0, 0) = 0.2;
  coefficients.alpha(0, 1) = 0.9;
  coefficients.alpha(1, 0) = 0.8;
  coefficients.alpha(1, 1) = 0.1;
  coefficients.beta(0, 0) = 0.2;
  coefficients.beta(0, 1) = 0.9;
  coefficients.beta(1, 0) = 0.8;
  coefficients.beta(1, 1) = 0.1;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  EXPECT_FALSE(nTest->notified);
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> controller(orbset, unrestrictedOccupations);
  controller.addSensitiveObject(nTest);
  auto dmat = controller.getDensityMatrix();
  EXPECT_DOUBLE_EQ(0.04, dmat.alpha(0, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.alpha(1, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.alpha(0, 1));
  EXPECT_DOUBLE_EQ(0.64, dmat.alpha(1, 1));
  EXPECT_DOUBLE_EQ(0.04, dmat.beta(0, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.beta(1, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.beta(0, 1));
  EXPECT_DOUBLE_EQ(0.64, dmat.beta(1, 1));
  controller.updateDensityMatrix();
  EXPECT_FALSE(nTest->notified);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Attach an OrbitalController (restricted)
 */
TEST_F(DensityMatrixControllerTest, AttachOrbitalController_OccVec_U) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::UNRESTRICTED>>();
  EXPECT_FALSE(nTest->notified);
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> controller(tRestrictedDensityMatrix);
  controller.addSensitiveObject(nTest);

  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficients = orbset->getCoefficients();
  coefficients.alpha(0, 0) = 0.1;
  coefficients.alpha(1, 0) = 0.9;
  coefficients.alpha(0, 1) = 0.8;
  coefficients.alpha(1, 1) = 0.2;
  coefficients.beta(0, 0) = 0.1;
  coefficients.beta(1, 0) = 0.9;
  coefficients.beta(0, 1) = 0.8;
  coefficients.beta(1, 1) = 0.2;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());

  controller.attachOrbitals(orbset, unrestrictedOccupations, true);

  auto dmat = controller.getDensityMatrix();
  EXPECT_FALSE(nTest->notified);
  EXPECT_DOUBLE_EQ(0.01, dmat.alpha(0, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.alpha(1, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.alpha(0, 1));
  EXPECT_DOUBLE_EQ(0.81, dmat.alpha(1, 1));
  EXPECT_DOUBLE_EQ(0.01, dmat.beta(0, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.beta(1, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.beta(0, 1));
  EXPECT_DOUBLE_EQ(0.81, dmat.beta(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Update via setter (unrestricted)
 */
TEST_F(DensityMatrixControllerTest, MatrixSetter_U) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::UNRESTRICTED>>();
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> controller(tUnrestrictedDensityMatrix);
  controller.addSensitiveObject(nTest);

  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> newMat(tBasisController);
  newMat.alpha(0, 0) = 1.01;
  newMat.alpha(0, 1) = 2.09;
  newMat.alpha(1, 0) = 2.09;
  newMat.alpha(1, 1) = 4.18;
  newMat.beta(0, 0) = 1.01;
  newMat.beta(0, 1) = 2.09;
  newMat.beta(1, 0) = 2.09;
  newMat.beta(1, 1) = 4.18;

  EXPECT_FALSE(nTest->notified);
  controller.setDensityMatrix(newMat);
  EXPECT_TRUE(nTest->notified);
  auto dmat = controller.getDensityMatrix();
  EXPECT_DOUBLE_EQ(1.01, dmat.alpha(0, 0));
  EXPECT_DOUBLE_EQ(2.09, dmat.alpha(1, 0));
  EXPECT_DOUBLE_EQ(2.09, dmat.alpha(0, 1));
  EXPECT_DOUBLE_EQ(4.18, dmat.alpha(1, 1));
  EXPECT_DOUBLE_EQ(1.01, dmat.beta(0, 0));
  EXPECT_DOUBLE_EQ(2.09, dmat.beta(1, 0));
  EXPECT_DOUBLE_EQ(2.09, dmat.beta(0, 1));
  EXPECT_DOUBLE_EQ(4.18, dmat.beta(1, 1));
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Update via notification (unrestricted)
 */
TEST_F(DensityMatrixControllerTest, NotifiedUpdate_U) {
  auto nTest = std::make_shared<DensityMatrixControllerNotificationTest<Options::SCF_MODES::UNRESTRICTED>>();
  EXPECT_FALSE(nTest->notified);
  auto orbset = std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(tBasisController, 0);
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficients = orbset->getCoefficients();
  // set initial coefficients
  coefficients.alpha(0, 0) = 0.2;
  coefficients.alpha(0, 1) = 0.9;
  coefficients.alpha(1, 0) = 0.8;
  coefficients.alpha(1, 1) = 0.1;
  coefficients.beta(0, 0) = 0.2;
  coefficients.beta(0, 1) = 0.9;
  coefficients.beta(1, 0) = 0.8;
  coefficients.beta(1, 1) = 0.1;
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  DensityMatrixController<Options::SCF_MODES::UNRESTRICTED> controller(orbset, unrestrictedOccupations);
  controller.addSensitiveObject(nTest);

  // check initial
  auto dmat = controller.getDensityMatrix();
  EXPECT_DOUBLE_EQ(0.04, dmat.alpha(0, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.alpha(1, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.alpha(0, 1));
  EXPECT_DOUBLE_EQ(0.64, dmat.alpha(1, 1));
  EXPECT_DOUBLE_EQ(0.04, dmat.beta(0, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.beta(1, 0));
  EXPECT_DOUBLE_EQ(0.16, dmat.beta(0, 1));
  EXPECT_DOUBLE_EQ(0.64, dmat.beta(1, 1));

  // edit obset
  coefficients.alpha(0, 0) = 0.1;
  coefficients.alpha(1, 0) = 0.9;
  coefficients.alpha(0, 1) = 0.8;
  coefficients.alpha(1, 1) = 0.2;
  coefficients.beta(0, 0) = 0.1;
  coefficients.beta(1, 0) = 0.9;
  coefficients.beta(0, 1) = 0.8;
  coefficients.beta(1, 1) = 0.2;
  EXPECT_FALSE(nTest->notified);

  // trigger orbset notification by setting new coefficients
  orbset->updateOrbitals(coefficients, orbset->getEigenvalues());
  dmat = controller.getDensityMatrix();
  EXPECT_TRUE(nTest->notified);
  EXPECT_DOUBLE_EQ(0.01, dmat.alpha(0, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.alpha(1, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.alpha(0, 1));
  EXPECT_DOUBLE_EQ(0.81, dmat.alpha(1, 1));
  EXPECT_DOUBLE_EQ(0.01, dmat.beta(0, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.beta(1, 0));
  EXPECT_DOUBLE_EQ(0.09, dmat.beta(0, 1));
  EXPECT_DOUBLE_EQ(0.81, dmat.beta(1, 1));
}
} /* end namespace Serenity */
