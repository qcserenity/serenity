/**
 * @file Ao2MoTransformer_test.cpp
 *
 * @date Nov 7, 2016
 * @author David Schnieders
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
#include "integrals/transformer/Ao2MoTransformer.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class Ao2MoTransformerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
/**
 * @test
 * @brief Analytical calculation of the multipole moments of H2 based on the density on a grid
 *
 * unrestricted, minimal basis
 */
TEST_F(Ao2MoTransformerTest, H2_Unrestricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto basisController = system->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, 1E-10);
  RegularRankFourTensor<double> eris(nBasisFunc);
  Ao2MoTransformer aoToMo(basisController);
  RegularRankFourTensor<double> result(nBasisFunc, 0.0);
  RegularRankFourTensor<double> resultSmall(3, 0.0);
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

  aoToMo.transformTwoElectronIntegrals(eris, result, coeffs, nBasisFunc);
  looper.loop(storeERIS);
  aoToMo.transformTwoElectronIntegrals(eris, resultSmall, coeffs, 3);
  /*
   * Check results
   */
  EXPECT_NEAR(-0.005698247, result(0, 7, 5, 0), 1e-5);
  EXPECT_NEAR(-0.011571752, result(1, 2, 6, 6), 1e-5);
  EXPECT_EQ(resultSmall(1, 2, 2, 2), result(1, 2, 2, 2));
  EXPECT_EQ(resultSmall(1, 1, 2, 1), result(1, 1, 2, 1));
};

} // namespace Serenity
