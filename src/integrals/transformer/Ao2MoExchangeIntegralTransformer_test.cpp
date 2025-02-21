/**
 * @file   Ao2MoExchangeIntegralTransformer_test.cpp
 * @author Moritz Bensberg
 * @date   25. December 2018
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
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //To be tested.
#include "analysis/PAOSelection/PNOConstructor.h"                   //Virtual basis.
#include "basis/BasisController.h"                                  //References.
#include "data/ElectronicStructure.h"                               //Access to orbital coefficients.
#include "data/OrbitalController.h"                                 //Access to orbital coefficients.
#include "data/OrbitalPair.h"                                       //Integrals are transformed pair-wise.
#include "data/OrbitalPairSet.h"                                    //OrbitalPairSet definition.
#include "data/PAOController.h"                                     //AO-PNO coefficients for reference transformation.
#include "integrals/looper/TwoElecFourCenterIntLooper.h"            //References.
#include "integrals/transformer/Ao2MoTransformer.h"                 //References.
#include "math/RegularRankFourTensor.h"                             //References.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"             //Access to k-set integrals.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                   //Access to kl-set integrals.
#include "postHF/LocalCorrelation/LocalCorrelationController.h"     //Easy construction of orbital pairs etc.
#include "system/SystemController.h"                                //Test systems.
#include "tasks/LocalizationTask.h"                                 //Orbital localization.
#include "tasks/ScfTask.h"                                          //Clean electronic structures.
#include "testsupply/SystemController__TEST_SUPPLY.h"               //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {
/**
 * @class Ao2MoExchangeIntegralTransformerTest Ao2MoExchangeIntegralTransformer_test.cpp
 * @brief Sets up the test of the Ao2MoExchangeIntegralTransformer.
 *
 * Note: The class Ao2MoExchangeIntegralTransformer is indirectly tested by the tests in LocalMP2Test,
 *       DLPNO_CCSD etc.
 */
class Ao2MoExchangeIntegralTransformerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }

  void calculate4CenterMOInts(std::shared_ptr<BasisController> aoBasis, const Eigen::MatrixXd& coeffs,
                              RegularRankFourTensor<double>& result) {
    // Calculte fully transformed integrals
    const unsigned int nBasisFunc = aoBasis->getNBasisFunctions();
    TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, aoBasis, 1E-10);
    RegularRankFourTensor<double> eris(nBasisFunc);
    Ao2MoTransformer aoToMo(aoBasis);
    auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                   const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
      (void)threadId; // no warning, please
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
    aoToMo.transformTwoElectronIntegrals(eris, result, coeffs, coeffs.cols());
  }
};

TEST_F(Ao2MoExchangeIntegralTransformerTest, H2_RIvsAO_allIntegrals) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  // Build reference
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_DEF2_TZVP);
  auto basisController = system->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();

  // Calculate exchange integrals (ia|jb)
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  auto pnoConstructor = localCorrelationController->producePNOConstructor();
  localCorrelationController->initializeSingles();
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs =
      localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
      system->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
      localCorrelationController->getMO3CenterIntegralController(), orbitalPairs, pnoConstructor);
  localCorrelationController->selectDistantOrbitalPairs();
  localCorrelationController->buildOrbitalPairCouplingMap();
  localCorrelationController->buildKLOrbitalPairs();

  localCorrelationController->getCloseOrbitalPairSets();
  Ao2MoExchangeIntegralTransformer::transformAllIntegrals(system->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
                                                          localCorrelationController->getMO3CenterIntegralController(),
                                                          orbitalPairs, false);

  auto pair = orbitalPairs[0];
  Eigen::MatrixXd R_00 = localCorrelationController->getPAOController()->getAllPAOs();

  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const Eigen::MatrixXd pnoCoefficients = R_00 * pair->toPAODomain;
  coeffs.block(0, 1, 12, 11) = pnoCoefficients;
  RegularRankFourTensor<double> result(nBasisFunc, 0.0);
  calculate4CenterMOInts(basisController, coeffs, result);
  unsigned int nPNOs = pair->ac_bd->rows();

  //(ia|jb)
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b < nPNOs; ++b) {
      double diff = pair->k_ij(a, b) - result(a + 1, 0, b + 1, 0);
      EXPECT_NEAR(0.0, diff, 5e-4);
    } // for b
  }   // for a

  //(ij|ab)
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b < nPNOs; ++b) {
      double diff = pair->ij_ab(a, b) - result(0, 0, a + 1, b + 1);
      EXPECT_NEAR(0.0, diff, 5e-4);
    } // for b
  }   // for a

  //(ib|ac),(jb|ac)
  // They should all be identical because i=j=0 and therefore [ij] = [i] = [j]
  for (unsigned int a = 0; a < nPNOs; ++a) {
    double symCheck = (pair->ia_bc[a] - pair->ja_bc[a]).array().abs().sum();
    EXPECT_NEAR(0.0, symCheck, 1e-9);
    for (unsigned int b = 0; b < nPNOs; ++b) {
      for (unsigned int c = 0; c < nPNOs; ++c) {
        double diff = pair->ia_bc[c](a, b) - result(0, a + 1, b + 1, c + 1);
        EXPECT_NEAR(0.0, diff, 8e-4);
      } // for c
    }   // for b
  }     // for a
  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
};

} /* namespace Serenity */
