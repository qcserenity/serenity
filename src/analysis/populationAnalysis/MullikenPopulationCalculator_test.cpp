/**
 * @file   MullikenPopulationCalculator_test.cpp
 *
 * @date   Sep 16, 2014
 * @author Jan Unsleber
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
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "math/Matrix.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory>

namespace Serenity {

class MullikenPopulationCalculatorTest : public ::testing::Test {
 protected:
};
/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(MullikenPopulationCalculatorTest, atomPopulationsUnrestricted) {
  /*
   * Set up test environment
   */
  auto basisController = BasisController__TEST_SUPPLY::getAtomCenteredBasisController(TEST_BASIS_CONTROLLERS::MINIMAL);
  /*
   * Set up input data
   */
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> densityMatrix(basisController);
  densityMatrix.alpha(0, 0) = 20.0;
  densityMatrix.alpha(0, 1) = 50.0;
  densityMatrix.alpha(1, 0) = 50.0;
  densityMatrix.alpha(1, 1) = 60.0;
  densityMatrix.beta(0, 0) = 20.0;
  densityMatrix.beta(0, 1) = 50.0;
  densityMatrix.beta(1, 0) = 50.0;
  densityMatrix.beta(1, 1) = 60.0;
  MatrixInBasis<RESTRICTED> overlaps(basisController);
  overlaps(0, 0) = 2.0;
  overlaps(0, 1) = 2.0;
  overlaps(1, 0) = 2.0;
  overlaps(1, 1) = 2.0;
  const auto indices = basisController->getBasisIndices();
  /*
   * Calculate results
   */
  MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED> Calculator;
  auto populations = Calculator.calculateAtomPopulations((const DensityMatrix<Options::SCF_MODES::UNRESTRICTED>)densityMatrix,
                                                         (const MatrixInBasis<RESTRICTED>)overlaps, indices);
  /*
   * Check results
   */
  EXPECT_EQ(populations.alpha[0], 140.0);
  EXPECT_EQ(populations.alpha[1], 220.0);
  EXPECT_EQ(populations.beta[0], 140.0);
  EXPECT_EQ(populations.beta[1], 220.0);
  EXPECT_NE(populations.alpha[0], 1.0);
  EXPECT_NE(populations.alpha[1], 1.0);
}
/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(MullikenPopulationCalculatorTest, atomPopulationsRestricted) {
  /*
   * Set up test environment
   */
  auto basisController = BasisController__TEST_SUPPLY::getAtomCenteredBasisController(TEST_BASIS_CONTROLLERS::MINIMAL);
  /*
   * Set up input data
   */
  DensityMatrix<Options::SCF_MODES::RESTRICTED> densityMatrix(basisController);
  densityMatrix(0, 0) = 20.0;
  densityMatrix(0, 1) = 50.0;
  densityMatrix(1, 0) = 50.0;
  densityMatrix(1, 1) = 60.0;
  MatrixInBasis<RESTRICTED> overlaps(basisController);
  overlaps(0, 0) = 2.0;
  overlaps(0, 1) = 2.0;
  overlaps(1, 0) = 2.0;
  overlaps(1, 1) = 2.0;
  const auto indices = basisController->getBasisIndices();
  /*
   * Calculate results
   */
  MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED> Calculator;
  auto populations = Calculator.calculateAtomPopulations((const DensityMatrix<Options::SCF_MODES::RESTRICTED>)densityMatrix,
                                                         (const MatrixInBasis<RESTRICTED>)overlaps, indices);
  /*
   * Check results
   */
  EXPECT_EQ(populations[0], 140.0);
  EXPECT_EQ(populations[1], 220.0);
  EXPECT_NE(populations[0], 1.0);
  EXPECT_NE(populations[1], 1.0);
}
/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(MullikenPopulationCalculatorTest, atomOrbPopulationsRestricted) {
  /*
   * Set up test environment
   */
  auto basisController = BasisController__TEST_SUPPLY::getAtomCenteredBasisController(TEST_BASIS_CONTROLLERS::MINIMAL);
  /*
   * Set up input data
   */
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffMatrix(basisController);
  coeffMatrix(0, 0) = 0.5;
  coeffMatrix(0, 1) = 0.5;
  coeffMatrix(1, 0) = 0.5;
  coeffMatrix(1, 1) = -0.5;
  MatrixInBasis<RESTRICTED> overlaps(basisController);
  overlaps(0, 0) = 1.0;
  overlaps(0, 1) = 0.1;
  overlaps(1, 0) = 0.1;
  overlaps(1, 1) = 1.0;
  const auto indices = basisController->getBasisIndices();
  /*
   * Calculate results
   */
  MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED> Calculator;
  auto populations =
      Calculator.calculateAtomwiseOrbitalPopulations((const CoefficientMatrix<Options::SCF_MODES::RESTRICTED>)coeffMatrix,
                                                     (const MatrixInBasis<RESTRICTED>)overlaps, indices);
  /*
   * Check results
   */
  EXPECT_EQ(populations(0, 0), 0.275);
  EXPECT_EQ(populations(0, 1), 0.225);
  EXPECT_EQ(populations(1, 0), 0.275);
  EXPECT_EQ(populations(1, 1), 0.225);
  EXPECT_NE(populations(0, 0), 1.0);
  EXPECT_NE(populations(0, 0), 1.0);
}
/**
 * @test
 * @brief Checks correct functionality for a tiny unphysical test case
 */
TEST_F(MullikenPopulationCalculatorTest, atomOrbPopulationsUnrestricted) {
  /*
   * Set up test environment
   */
  auto basisController = BasisController__TEST_SUPPLY::getAtomCenteredBasisController(TEST_BASIS_CONTROLLERS::MINIMAL);
  /*
   * Set up input data
   */
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coeffMatrix(basisController);
  coeffMatrix.alpha(0, 0) = 0.5;
  coeffMatrix.alpha(0, 1) = 0.5;
  coeffMatrix.alpha(1, 0) = 0.5;
  coeffMatrix.alpha(1, 1) = -0.5;
  coeffMatrix.beta(0, 0) = 0.5;
  coeffMatrix.beta(0, 1) = 0.5;
  coeffMatrix.beta(1, 0) = 0.5;
  coeffMatrix.beta(1, 1) = -0.5;
  MatrixInBasis<RESTRICTED> overlaps(basisController);
  overlaps(0, 0) = 1.0;
  overlaps(0, 1) = 0.1;
  overlaps(1, 0) = 0.1;
  overlaps(1, 1) = 1.0;
  const auto indices = basisController->getBasisIndices();
  /*
   * Calculate results
   */
  MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED> Calculator;
  auto populations =
      Calculator.calculateAtomwiseOrbitalPopulations((const CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>)coeffMatrix,
                                                     (const MatrixInBasis<RESTRICTED>)overlaps, indices);
  /*
   * Check results
   */
  EXPECT_EQ(populations.alpha(0, 0), 0.275);
  EXPECT_EQ(populations.alpha(0, 1), 0.225);
  EXPECT_EQ(populations.alpha(1, 0), 0.275);
  EXPECT_EQ(populations.alpha(1, 1), 0.225);
  EXPECT_EQ(populations.beta(0, 0), 0.275);
  EXPECT_EQ(populations.beta(0, 1), 0.225);
  EXPECT_EQ(populations.beta(1, 0), 0.275);
  EXPECT_EQ(populations.beta(1, 1), 0.225);
  EXPECT_NE(populations.beta(0, 0), 1.0);
  EXPECT_NE(populations.beta(0, 0), 1.0);
}
} /* namespace Serenity */
