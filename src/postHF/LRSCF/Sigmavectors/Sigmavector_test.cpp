/**
 * @file Sigmavector_test.cpp
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers
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

/* Include Class Header*/
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"

/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/transformer/Ao2MoTransformer.h"
#include "math/RegularRankFourTensor.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/FockSigmavector.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @class SigmavectorTest
 * @brief Sets everything up for the tests of Sigmavector.h/.cpp .
 */
class SigmavectorTest : public ::testing::Test {
 protected:
  SigmavectorTest() {
  }

  virtual ~SigmavectorTest() = default;

  /// system
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

TEST_F(SigmavectorTest, DJK_res) {
  // Setup test system
  auto systemController = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP_CAMB3LYP);
  systemController->getElectronicStructure<Options::SCF_MODES::RESTRICTED>();
  auto basisController = systemController->getBasisController();

  // CAMB3LYP settings
  double exchangeRatio = 0.19;
  double lrexchangeRatio = 0.46;
  double mu = 0.33;

  // Calculate AO integrals and transform to MO basis
  RegularRankFourTensor<double> eris(basisController->getNBasisFunctions(), 0.0);
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, 1E-10);
  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warnings, please
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };
  looper.loop(storeERIS);

  RegularRankFourTensor<double> lreris(basisController->getNBasisFunctions(), 0.0);
  TwoElecFourCenterIntLooper lrlooper(LIBINT_OPERATOR::erf_coulomb, 0, basisController, 1E-10, mu);
  auto const storeLRERIS = [&lreris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                     const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warnings, please
    lreris(b, a, i, j) = integral(0);
    lreris(b, a, j, i) = integral(0);
    lreris(a, b, j, i) = integral(0);
    lreris(a, b, i, j) = integral(0);
    lreris(i, j, b, a) = integral(0);
    lreris(i, j, a, b) = integral(0);
    lreris(j, i, b, a) = integral(0);
    lreris(j, i, a, b) = integral(0);
  };
  lrlooper.loop(storeLRERIS);

  Ao2MoTransformer ao2mo(basisController);
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients =
      systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  ao2mo.transformTwoElectronIntegrals(eris, eris, coefficients, basisController->getNBasisFunctions());
  ao2mo.transformTwoElectronIntegrals(lreris, lreris, coefficients, basisController->getNBasisFunctions());

  // In order to test the subsystem sigma vectors we need intersubsystem four center MOs. Since the
  // transformation for exchange and Coulomb integrals differs, this gets messy. We can however test
  // the subsytem sigma vectors by using the same test system twice. Physically, this is nonsense but
  // this trick can be justified since all parts of the code are tested.
  auto nOcc = systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  auto nVirt = systemController->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nDimension = nOcc * nVirt;
  const auto orbitalEnergies =
      systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();

  Eigen::MatrixXd J(nDimension, nDimension);
  J.setZero();

  Eigen::MatrixXd D(nDimension, nDimension);
  D.setZero();

  Eigen::MatrixXd KP(nDimension, nDimension);
  KP.setZero();

  for (unsigned int i = 0, ia = 0; i < nOcc; ++i) {
    for (unsigned int a = nOcc; a < nOcc + nVirt; ++a, ++ia) {
      for (unsigned int j = 0, jb = 0; j < nOcc; ++j) {
        for (unsigned int b = nOcc; b < nOcc + nVirt; ++b, ++jb) {
          J(ia, jb) += eris(a, i, b, j);
          if (i == j && a == b)
            D(ia, jb) += orbitalEnergies(a) - orbitalEnergies(i);
          KP(ia, jb) += exchangeRatio * (eris(a, j, b, i) + eris(a, b, j, i));
          KP(ia, jb) += lrexchangeRatio * (lreris(a, j, b, i) + lreris(a, b, j, i));
        }
      }
    }
  }

  // Setup two sets of guess vector
  std::vector<Eigen::MatrixXd> guess(2);
  guess[0].resize(2 * nDimension, 2);
  guess[1].resize(2 * nDimension, 2);
  guess[0].setRandom();
  guess[1].setRandom();

  // Setup reference orbitals
  LRSCFTaskSettings settings;
  settings.embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::CAMB3LYP;
  std::shared_ptr<LRSCFController<Options::SCF_MODES::RESTRICTED>> lrscfController =
      std::make_shared<LRSCFController<Options::SCF_MODES::RESTRICTED>>(systemController, settings);

  // Calculate Sigma vectors
  std::vector<Eigen::MatrixXd> J_correct(guess.size());

  std::vector<Eigen::MatrixXd> D_correct(guess.size());

  std::vector<Eigen::MatrixXd> KP_correct(guess.size());

  for (unsigned int iSet = 0; iSet < guess.size(); ++iSet) {
    J_correct[iSet].resize(2 * nDimension, guess[iSet].cols());
    J_correct[iSet].setZero();
    // sigma_I = sum_J M_IJ*b_J (M_IJ is the same for all I, J since subsystem I=J)
    J_correct[iSet].block(0, 0, nDimension, guess[iSet].cols()) +=
        J * guess[iSet].block(0, 0, nDimension, guess[iSet].cols());
    J_correct[iSet].block(0, 0, nDimension, guess[iSet].cols()) +=
        J * guess[iSet].block(nDimension, 0, nDimension, guess[iSet].cols());
    // Sigma vectors of both subsystems are equal...
    J_correct[iSet].block(nDimension, 0, nDimension, guess[iSet].cols()) =
        J_correct[iSet].block(0, 0, nDimension, guess[iSet].cols());

    D_correct[iSet].resize(2 * nDimension, guess[iSet].cols());
    D_correct[iSet].setZero();
    D_correct[iSet].block(0, 0, nDimension, guess[iSet].cols()) +=
        D * guess[iSet].block(0, 0, nDimension, guess[iSet].cols());
    D_correct[iSet].block(nDimension, 0, nDimension, guess[iSet].cols()) =
        D * guess[iSet].block(nDimension, 0, nDimension, guess[iSet].cols());

    KP_correct[iSet].resize(2 * nDimension, guess[iSet].cols());
    KP_correct[iSet].setZero();
    KP_correct[iSet].block(0, 0, nDimension, guess[iSet].cols()) +=
        KP * guess[iSet].block(0, 0, nDimension, guess[iSet].cols());
    KP_correct[iSet].block(0, 0, nDimension, guess[iSet].cols()) +=
        KP * guess[iSet].block(nDimension, 0, nDimension, guess[iSet].cols());
    KP_correct[iSet].block(nDimension, 0, nDimension, guess[iSet].cols()) =
        KP_correct[iSet].block(0, 0, nDimension, guess[iSet].cols());
  }

  // Test Coulomb contribution
  CoulombSigmavector<Options::SCF_MODES::RESTRICTED> JSigmavector({lrscfController, lrscfController}, guess);
  std::vector<Eigen::MatrixXd> J_test = JSigmavector.getSigma();
  // Test orbital energy difference contribution
  FockSigmavector<Options::SCF_MODES::RESTRICTED> DSigmavector({lrscfController, lrscfController}, guess);
  std::vector<Eigen::MatrixXd> D_test = DSigmavector.getSigma();
  // Test exchange contributions
  ExchangeSigmavector<Options::SCF_MODES::RESTRICTED> KPSigmavector({lrscfController, lrscfController}, guess, {1, 1},
                                                                    false, false);
  std::vector<Eigen::MatrixXd> KP_test = KPSigmavector.getSigma();

  for (unsigned int iSet = 0; iSet < guess.size(); ++iSet) {
    for (unsigned int iGuess = 0; iGuess < guess[iSet].cols(); ++iGuess) {
      for (unsigned int ia = 0; ia < 2 * nDimension; ++ia) {
        EXPECT_NEAR(J_correct[iSet](ia, iGuess), J_test[iSet](ia, iGuess), 1.0e-7);
        EXPECT_NEAR(D_correct[iSet](ia, iGuess), D_test[iSet](ia, iGuess), 1.0e-7);
        EXPECT_NEAR(KP_correct[iSet](ia, iGuess), KP_test[iSet](ia, iGuess), 1.0e-7);
      }
    }
  }
}

} // namespace Serenity
