/**
 * @file   Damper_test.cpp
 *
 * @date   Feb 07, 2024
 * @author Lukas Paetow
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
#include "scf/damper/Damper.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class DamperTest : public ::testing::Test {
 protected:
  DamperTest() : tBasisController(BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::MINIMAL)) {
  }
  virtual ~DamperTest() = default;
  std::shared_ptr<BasisController> tBasisController;
};

/**
 * @test
 * @brief Tests Damper.h: Restricted Static
 */
TEST_F(DamperTest, Restricted_StaticDamping) {
  FockMatrix<RESTRICTED> matrix(tBasisController);
  matrix(0, 0) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(0, 1) = 3.0;
  matrix(1, 1) = 4.0;

  FockMatrix<RESTRICTED> matrix2(tBasisController);
  matrix2(0, 0) = 9.0;
  matrix2(1, 0) = 8.0;
  matrix2(0, 1) = 7.0;
  matrix2(1, 1) = 6.0;

  Damper<Options::SCF_MODES::RESTRICTED> damper(0.5);
  damper.staticDamp(matrix);
  damper.staticDamp(matrix2);

  EXPECT_EQ(matrix2(0, 0), 5.0);
  EXPECT_EQ(matrix2(1, 0), 5.0);
  EXPECT_EQ(matrix2(0, 1), 5.0);
  EXPECT_EQ(matrix2(1, 1), 5.0);
}

/**
 * @test
 * @brief Tests Damper.h: Restricted Dynamic
 */
TEST_F(DamperTest, Restricted_DynamicDamping) {
  FockMatrix<RESTRICTED> matrix(tBasisController);
  matrix(0, 0) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(0, 1) = 3.0;
  matrix(1, 1) = 4.0;

  DensityMatrix<RESTRICTED> matrix2(tBasisController);
  matrix2(0, 0) = 9.0;
  matrix2(1, 0) = 8.0;
  matrix2(0, 1) = 7.0;
  matrix2(1, 1) = 6.0;

  Damper<Options::SCF_MODES::RESTRICTED> damper;
  damper.dynamicDamp(matrix, matrix2);

  matrix(0, 0) = 2.0;
  matrix(1, 0) = 3.0;
  matrix(0, 1) = 4.0;
  matrix(1, 1) = 5.0;
  matrix2(0, 0) = 10.0;
  matrix2(1, 0) = 9.0;
  matrix2(0, 1) = 8.0;
  matrix2(1, 1) = 7.0;

  damper.dynamicDamp(matrix, matrix2);

  EXPECT_EQ(matrix(0, 0), -0.25);
  EXPECT_EQ(matrix(1, 0), 0.75);
  EXPECT_EQ(matrix(0, 1), 1.75);
  EXPECT_EQ(matrix(1, 1), 2.75);
}

/**
 * @test
 * @brief Tests Damper.h: Restricted Arithmetic Series
 */
TEST_F(DamperTest, Restricted_ArithmeticSeriesDamping) {
  FockMatrix<RESTRICTED> matrix(tBasisController);
  matrix(0, 0) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(0, 1) = 3.0;
  matrix(1, 1) = 4.0;

  FockMatrix<RESTRICTED> matrix2(tBasisController);
  matrix2(0, 0) = 9.0;
  matrix2(1, 0) = 8.0;
  matrix2(0, 1) = 7.0;
  matrix2(1, 1) = 6.0;

  FockMatrix<RESTRICTED> matrix3(tBasisController);
  matrix3(0, 0) = 10.0;
  matrix3(1, 0) = 9.0;
  matrix3(0, 1) = 8.0;
  matrix3(1, 1) = 7.0;

  double dStart = 0.7;
  double dStep = 0.05;
  double dEnd = 0.2;
  int iStartUp = 2;
  Damper<Options::SCF_MODES::RESTRICTED> damper(dStart, dStep, dEnd, iStartUp);
  damper.arithmeticSeriesDamp(matrix);
  damper.arithmeticSeriesDamp(matrix2);
  damper.arithmeticSeriesDamp(matrix3);
  damper.arithmeticSeriesDamp(matrix2);

  EXPECT_NEAR(matrix2(0, 0), 4.687, 1e-8);
  EXPECT_NEAR(matrix2(1, 0), 4.814, 1e-8);
  EXPECT_NEAR(matrix2(0, 1), 4.941, 1e-8);
  EXPECT_NEAR(matrix2(1, 1), 5.068, 1e-8);
}

/**
 * @test
 * @brief Tests Damper.h: Unrestricted Static
 */
TEST_F(DamperTest, Unrestricted_StaticDamping) {
  FockMatrix<UNRESTRICTED> matrix(tBasisController);
  matrix.alpha(0, 0) = 1.0;
  matrix.alpha(1, 0) = 2.0;
  matrix.alpha(0, 1) = 3.0;
  matrix.alpha(1, 1) = 4.0;
  matrix.beta(0, 0) = 1.0;
  matrix.beta(1, 0) = 2.0;
  matrix.beta(0, 1) = 3.0;
  matrix.beta(1, 1) = 4.0;

  FockMatrix<UNRESTRICTED> matrix2(tBasisController);
  matrix2.alpha(0, 0) = 9.0;
  matrix2.alpha(1, 0) = 8.0;
  matrix2.alpha(0, 1) = 7.0;
  matrix2.alpha(1, 1) = 6.0;
  matrix2.beta(0, 0) = 9.0;
  matrix2.beta(1, 0) = 8.0;
  matrix2.beta(0, 1) = 7.0;
  matrix2.beta(1, 1) = 6.0;

  Damper<Options::SCF_MODES::UNRESTRICTED> damper(0.5);
  damper.staticDamp(matrix);
  damper.staticDamp(matrix2);

  EXPECT_EQ(matrix2.alpha(0, 0), 5.0);
  EXPECT_EQ(matrix2.alpha(1, 0), 5.0);
  EXPECT_EQ(matrix2.alpha(0, 1), 5.0);
  EXPECT_EQ(matrix2.alpha(1, 1), 5.0);

  EXPECT_EQ(matrix2.beta(0, 0), 5.0);
  EXPECT_EQ(matrix2.beta(1, 0), 5.0);
  EXPECT_EQ(matrix2.beta(0, 1), 5.0);
  EXPECT_EQ(matrix2.beta(1, 1), 5.0);
}

/**
 * @test
 * @brief Tests Damper.h: Unrestricted Dynamic
 */
TEST_F(DamperTest, Unrestricted_DynamicDamping) {
  FockMatrix<UNRESTRICTED> matrix(tBasisController);
  matrix.alpha(0, 0) = 1.0;
  matrix.alpha(1, 0) = 2.0;
  matrix.alpha(0, 1) = 3.0;
  matrix.alpha(1, 1) = 4.0;
  matrix.beta(0, 0) = 1.0;
  matrix.beta(1, 0) = 2.0;
  matrix.beta(0, 1) = 3.0;
  matrix.beta(1, 1) = 4.0;

  DensityMatrix<UNRESTRICTED> matrix2(tBasisController);
  matrix2.alpha(0, 0) = 9.0;
  matrix2.alpha(1, 0) = 8.0;
  matrix2.alpha(0, 1) = 7.0;
  matrix2.alpha(1, 1) = 6.0;
  matrix2.beta(0, 0) = 9.0;
  matrix2.beta(1, 0) = 8.0;
  matrix2.beta(0, 1) = 7.0;
  matrix2.beta(1, 1) = 6.0;

  Damper<Options::SCF_MODES::UNRESTRICTED> damper;
  damper.dynamicDamp(matrix, matrix2);

  matrix.alpha(0, 0) = 2.0;
  matrix.alpha(1, 0) = 3.0;
  matrix.alpha(0, 1) = 4.0;
  matrix.alpha(1, 1) = 5.0;
  matrix.beta(0, 0) = 2.0;
  matrix.beta(1, 0) = 3.0;
  matrix.beta(0, 1) = 4.0;
  matrix.beta(1, 1) = 5.0;
  matrix2.alpha(0, 0) = 10.0;
  matrix2.alpha(1, 0) = 9.0;
  matrix2.alpha(0, 1) = 8.0;
  matrix2.alpha(1, 1) = 7.0;
  matrix2.beta(0, 0) = 10.0;
  matrix2.beta(1, 0) = 9.0;
  matrix2.beta(0, 1) = 8.0;
  matrix2.beta(1, 1) = 7.0;

  damper.dynamicDamp(matrix, matrix2);

  EXPECT_EQ(matrix.alpha(0, 0), -0.25);
  EXPECT_EQ(matrix.alpha(1, 0), 0.75);
  EXPECT_EQ(matrix.alpha(0, 1), 1.75);
  EXPECT_EQ(matrix.alpha(1, 1), 2.75);
  EXPECT_EQ(matrix.beta(0, 0), -0.25);
  EXPECT_EQ(matrix.beta(1, 0), 0.75);
  EXPECT_EQ(matrix.beta(0, 1), 1.75);
  EXPECT_EQ(matrix.beta(1, 1), 2.75);
}

/**
 * @test
 * @brief Tests Damper.h: Unrestricted Arithmetic Series
 */
TEST_F(DamperTest, Unrestricted_ArithmeticSeriesDamping) {
  FockMatrix<UNRESTRICTED> matrix(tBasisController);
  matrix.alpha(0, 0) = 1.0;
  matrix.alpha(1, 0) = 2.0;
  matrix.alpha(0, 1) = 3.0;
  matrix.alpha(1, 1) = 4.0;
  matrix.beta(0, 0) = 1.0;
  matrix.beta(1, 0) = 2.0;
  matrix.beta(0, 1) = 3.0;
  matrix.beta(1, 1) = 4.0;

  FockMatrix<UNRESTRICTED> matrix2(tBasisController);
  matrix2.alpha(0, 0) = 9.0;
  matrix2.alpha(1, 0) = 8.0;
  matrix2.alpha(0, 1) = 7.0;
  matrix2.alpha(1, 1) = 6.0;
  matrix2.beta(0, 0) = 9.0;
  matrix2.beta(1, 0) = 8.0;
  matrix2.beta(0, 1) = 7.0;
  matrix2.beta(1, 1) = 6.0;

  FockMatrix<UNRESTRICTED> matrix3(tBasisController);
  matrix3.alpha(0, 0) = 10.0;
  matrix3.alpha(1, 0) = 9.0;
  matrix3.alpha(0, 1) = 8.0;
  matrix3.alpha(1, 1) = 7.0;
  matrix3.beta(0, 0) = 10.0;
  matrix3.beta(1, 0) = 9.0;
  matrix3.beta(0, 1) = 8.0;
  matrix3.beta(1, 1) = 7.0;

  double dStart = 0.7;
  double dStep = 0.05;
  double dEnd = 0.2;
  int iStartUp = 2;

  Damper<Options::SCF_MODES::UNRESTRICTED> damper(dStart, dStep, dEnd, iStartUp);
  damper.arithmeticSeriesDamp(matrix);
  damper.arithmeticSeriesDamp(matrix2);
  damper.arithmeticSeriesDamp(matrix3);
  damper.arithmeticSeriesDamp(matrix2);

  EXPECT_NEAR(matrix2.alpha(0, 0), 4.687, 1e-8);
  EXPECT_NEAR(matrix2.alpha(1, 0), 4.814, 1e-8);
  EXPECT_NEAR(matrix2.alpha(0, 1), 4.941, 1e-8);
  EXPECT_NEAR(matrix2.alpha(1, 1), 5.068, 1e-8);

  EXPECT_NEAR(matrix2.beta(0, 0), 4.687, 1e-8);
  EXPECT_NEAR(matrix2.beta(1, 0), 4.814, 1e-8);
  EXPECT_NEAR(matrix2.beta(0, 1), 4.941, 1e-8);
  EXPECT_NEAR(matrix2.beta(1, 1), 5.068, 1e-8);
}

} /* namespace Serenity */
