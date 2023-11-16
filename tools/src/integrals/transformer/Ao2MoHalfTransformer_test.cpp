/**
 * @file   Ao2MoHalfTransformer_test.cpp
 * @author Moritz Bensberg
 * @date   17. December 2018
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
#include "integrals/transformer/Ao2MoHalfTransformer.h"  //To be tested.
#include "basis/BasisController.h"                       //Basis controller definition.
#include "data/ElectronicStructure.h"                    //Orbital coefficients.
#include "data/OrbitalController.h"                      //Orbital coefficients.
#include "integrals/looper/TwoElecFourCenterIntLooper.h" //References.
#include "integrals/transformer/Ao2MoTransformer.h"      //References.
#include "math/RegularRankFourTensor.h"                  //References.
#include "system/SystemController.h"                     //Test systems.
#include "testsupply/SystemController__TEST_SUPPLY.h"    //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class Ao2MoHalfTransformerTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
  @test
  @brief Tests the Ao2MoHalfTransformer by comparing to the fully transformed
         integrals to the integrals obtained after manual transformation of the
         half transformed integrals.
*/
TEST_F(Ao2MoHalfTransformerTest, CO_Restricted) {
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::CO_MINBAS);
  auto basisController = system->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  // Calculte fully transformed integrals
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, 1E-10);
  RegularRankFourTensor<double> eris(nBasisFunc);
  Ao2MoTransformer aoToMo(basisController);
  RegularRankFourTensor<double> result(nBasisFunc, 0.0);
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
  // Calculate half transformed integrals
  Ao2MoHalfTransformer halfTransformer(basisController, basisController);
  Eigen::MatrixXd halfTransformed;
  auto c_i = coeffs.col(0);
  Eigen::MatrixXd pairDensityMatrix = c_i * c_i.transpose();
  halfTransformer.transformTwoElectronIntegrals(halfTransformed, pairDensityMatrix);
  // Manual transformation
  Eigen::MatrixXd tmp = (coeffs.transpose() * halfTransformed * coeffs).eval();
  // Compare the integrals
  for (unsigned int row = 0; row < tmp.rows(); ++row) {
    for (unsigned int col = 0; col < tmp.cols(); ++col) {
      EXPECT_NEAR(result(row, 0, col, 0), tmp(row, col), 1e-9);
    }
  }
};

} // namespace Serenity
