/**
 * @file SigmaVector_test.cpp
 *
 * @date Oct 09, 2017
 * @author Michael Boeckers
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
#include "integrals/transformer/Ao2MoTransformer.h"
#include "postHF/LRSCF/SigmaVector/DeltaESigmaVector.h"
#include "postHF/LRSCF/SigmaVector/GSigmaVector.h"
#include "postHF/LRSCF/SigmaVector/KSigmaVector.h"
#include "data/OrbitalController.h"
#include "math/RegularRankFourTensor.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class SigmaVectorTest
 * @brief Sets everything up for the tests of SigmaVector.h/.cpp .
 */
class SigmaVectorTest : public ::testing::Test {
 protected:
  SigmaVectorTest() {
  }

  virtual ~SigmaVectorTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SigmaVectorTest,RTDHF) {
  //Setup test system
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto basisController = systemController->getBasisController();

  //Calculate AO integrals and transform to MO basis
  RegularRankFourTensor<double> eris(basisController->getNBasisFunctions(), 0.0);
  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb,0,basisController,1E-10);
  auto const storeERIS = [&eris]
                          (const unsigned int&  a,
                              const unsigned int&  b,
                              const unsigned int&  i,
                              const unsigned int&  j,
                              const Eigen::VectorXd&  integral,
                              const unsigned int threadId) {
    (void)threadId; //no warnings, please
    eris(b,a,i,j) = integral(0);
    eris(b,a,j,i) = integral(0);
    eris(a,b,j,i) = integral(0);
    eris(a,b,i,j) = integral(0);
    eris(i,j,b,a) = integral(0);
    eris(i,j,a,b) = integral(0);
    eris(j,i,b,a) = integral(0);
    eris(j,i,a,b) = integral(0);
  };
  looper.loop(storeERIS);
  Ao2MoTransformer ao2mo(basisController);
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients=systemController->
      getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  ao2mo.transformTwoElectronIntegrals(eris,eris,coefficients,basisController->getNBasisFunctions());
  //Setup A and B matrices (primitive)
  const auto orbitalEnergies = systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();
  auto nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto nVirt = systemController->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nDimension = nOcc*nVirt;
  Eigen::MatrixXd A(nDimension,nDimension);
  Eigen::MatrixXd B(nDimension,nDimension);
  A.setZero();
  B.setZero();
  for (unsigned int i = 0,ia = 0; i < nOcc; ++i) {
    for (unsigned int a = nOcc; a < nOcc + nVirt; ++a, ++ia) {
      for (unsigned int j = 0, jb = 0; j < nOcc; ++j) {
        for (unsigned int b = nOcc; b < nOcc + nVirt; ++b, ++jb) {
          if (i == j && a == b) A(ia,jb) += orbitalEnergies(a) - orbitalEnergies(i);
          A(ia,jb) += 2.0 * eris(a,i,b,j) - eris(a,b,j,i);
          B(ia,jb) += 2.0 * eris(i,a,j,b) - eris(a,j,b,i);
        }
      }
    }
  }

  //Setup guess vector
  Eigen::MatrixXd guess(nDimension,1);
  guess.setRandom();
  //Calculate differences between sigma vectors and check if zero
  DeltaESigmaVector<Options::SCF_MODES::RESTRICTED> deltaESigma(systemController,guess);
  KSigmaVector<Options::SCF_MODES::RESTRICTED> kSigma(systemController,guess,-1.0);
  GSigmaVector<Options::SCF_MODES::RESTRICTED> gSigma(systemController,guess);
  Eigen::VectorXd thisShouldBeZero = (A-B)*guess.col(0) - Eigen::MatrixXd(deltaESigma.getSigma() + kSigma.getSigma()).col(0);
  thisShouldBeZero += (A+B)*guess.col(0) - Eigen::MatrixXd(deltaESigma.getSigma() + gSigma.getSigma()).col(0);
  for (unsigned int i = 0; i < nDimension; ++i) {
    EXPECT_NEAR(0.0, thisShouldBeZero(i), 1.0e-7);
  }
}

}



