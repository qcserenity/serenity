/**
 * @file   CoefficientMatrix_test.cpp
 * @author Jan Unsleber
 *
 * @date   May 9, 2017
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
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "settings/Options.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
class CoefficientMatrixTest : public ::testing::Test {
 protected:
  CoefficientMatrixTest()
    : tBasisController(BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::MINIMAL)),
      tr(tBasisController),
      tu(tBasisController) {
    tr(0, 0) = 0.0;
    tr(0, 1) = 0.1;
    tr(1, 0) = 1.0;
    tr(1, 1) = 1.1;
    tu.alpha(0, 0) = 0.0;
    tu.alpha(0, 1) = 0.1;
    tu.alpha(1, 0) = 1.0;
    tu.alpha(1, 1) = 1.1;
    tu.beta(0, 0) = 0.0;
    tu.beta(0, 1) = 0.1;
    tu.beta(1, 0) = 1.0;
    tu.beta(1, 1) = 1.1;
  }
  virtual ~CoefficientMatrixTest() = default;
  // The basis in which it is defined
  std::shared_ptr<BasisController> tBasisController;
  // The object under test
  CoefficientMatrix<RESTRICTED> tr;
  CoefficientMatrix<UNRESTRICTED> tu;
};
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Normal construction.
 */
TEST_F(CoefficientMatrixTest, BasisConstruct_R) {
  CoefficientMatrix<RESTRICTED> matrix(tBasisController);
  // Check for zero initialization
  EXPECT_EQ(matrix(0, 0), 0.0);
  EXPECT_EQ(matrix(0, 1), 0.0);
  EXPECT_EQ(matrix(1, 0), 0.0);
  EXPECT_EQ(matrix(1, 1), 0.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Copy construction.
 */
TEST_F(CoefficientMatrixTest, CopyConstruct_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  // Check for zero initialization
  EXPECT_EQ(matrix(0, 0), 0.0);
  EXPECT_EQ(matrix(0, 1), 0.1);
  EXPECT_EQ(matrix(1, 0), 1.0);
  EXPECT_EQ(matrix(1, 1), 1.1);
  EXPECT_EQ(tr(0, 0), 0.0);
  EXPECT_EQ(tr(0, 1), 0.1);
  EXPECT_EQ(tr(1, 0), 1.0);
  EXPECT_EQ(tr(1, 1), 1.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= MatrixXd.
 */
TEST_F(CoefficientMatrixTest, Eq_MatrixXd_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  Eigen::MatrixXd test(2, 2);
  test << 1.0, 2.0, 2.0, 3.0;
  matrix = test;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.0);
  EXPECT_EQ(matrix(1, 0), 2.0);
  EXPECT_EQ(matrix(1, 1), 3.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= MiB.
 */
TEST_F(CoefficientMatrixTest, Eq_MatInB_3_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  MatrixInBasis<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  matrix = test;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.0);
  EXPECT_EQ(matrix(1, 0), 2.0);
  EXPECT_EQ(matrix(1, 1), 3.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, Eq_MatInB_1_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  CoefficientMatrix<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  matrix = test;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.0);
  EXPECT_EQ(matrix(1, 0), 2.0);
  EXPECT_EQ(matrix(1, 1), 3.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, Eq_Ref_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  matrix(0, 0) = 500.0;
  const auto& test = &matrix;
  EXPECT_EQ((*test)(0, 0), 500.0);
  EXPECT_EQ((*test)(0, 1), 0.1);
  EXPECT_EQ((*test)(1, 0), 1.0);
  EXPECT_EQ((*test)(1, 1), 1.1);
  EXPECT_EQ(test->getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, Eq_MatInB_2_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  auto bas2 = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  CoefficientMatrix<RESTRICTED> test(bas2);
  EXPECT_THROW(matrix.operator=(test), SerenityError);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator= Eigen3 computed expression.
 */
TEST_F(CoefficientMatrixTest, Eq_MatProd_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  Eigen::MatrixXd test1(2, 2);
  test1 << 1.0, 2.0, 2.0, 3.0;
  Eigen::MatrixXd test2(2, 2);
  test2 << 1.0, 0.0, 0.0, 1.0;
  matrix = test1 * test2;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.0);
  EXPECT_EQ(matrix(1, 0), 2.0);
  EXPECT_EQ(matrix(1, 1), 3.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator+= MatrixXd.
 */
TEST_F(CoefficientMatrixTest, PE_MatrixXd_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  Eigen::MatrixXd test(2, 2);
  test << 1.0, 2.0, 2.0, 3.0;
  matrix += test;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.1);
  EXPECT_EQ(matrix(1, 0), 3.0);
  EXPECT_EQ(matrix(1, 1), 4.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator-= MatrixXd.
 */
TEST_F(CoefficientMatrixTest, ME_MatrixXd_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  Eigen::MatrixXd test(2, 2);
  test << 1.0, 2.0, 2.0, 3.0;
  matrix -= test;
  EXPECT_EQ(matrix(0, 0), -1.0);
  EXPECT_EQ(matrix(0, 1), -1.9);
  EXPECT_EQ(matrix(1, 0), -1.0);
  EXPECT_EQ(matrix(1, 1), -1.9);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator+= CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, PE_MatInB_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  CoefficientMatrix<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  matrix += test;
  EXPECT_EQ(matrix(0, 0), 1.0);
  EXPECT_EQ(matrix(0, 1), 2.1);
  EXPECT_EQ(matrix(1, 0), 3.0);
  EXPECT_EQ(matrix(1, 1), 4.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator-= CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, ME_MatInB_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  CoefficientMatrix<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  matrix -= test;
  EXPECT_EQ(matrix(0, 0), -1.0);
  EXPECT_EQ(matrix(0, 1), -1.9);
  EXPECT_EQ(matrix(1, 0), -1.0);
  EXPECT_EQ(matrix(1, 1), -1.9);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator* CoefficientMatrix.
 */
TEST_F(CoefficientMatrixTest, Multiply_MatInB_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  CoefficientMatrix<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  auto result = matrix * test;
  EXPECT_EQ(result.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: operator* mixed.
 */
TEST_F(CoefficientMatrixTest, Multiply_Mixed_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  CoefficientMatrix<RESTRICTED> test1(tBasisController);
  test1(0, 0) = 1.0;
  test1(1, 0) = 2.0;
  test1(0, 1) = 2.0;
  test1(1, 1) = 3.0;
  Eigen::MatrixXd test2(2, 2);
  test2 << 1.0, 0.0, 0.0, 1.0;
  Eigen::MatrixXd result = test1 * test2;
  EXPECT_EQ(result(0, 0), 1.0);
  EXPECT_EQ(result(0, 1), 2.0);
  EXPECT_EQ(result(1, 0), 2.0);
  EXPECT_EQ(result(1, 1), 3.0);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: getBasisController(), getNBasisFunctions()
 */
TEST_F(CoefficientMatrixTest, GetBasis_R) {
  CoefficientMatrix<RESTRICTED> matrix(tBasisController);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
  EXPECT_EQ(matrix.getNBasisFunctions(), (unsigned int)2);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: total()
 */
TEST_F(CoefficientMatrixTest, Total_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  auto total = matrix.total();
  EXPECT_EQ(total(0, 0), 0.0);
  EXPECT_EQ(total(0, 1), 0.1);
  EXPECT_EQ(total(1, 0), 1.0);
  EXPECT_EQ(total(1, 1), 1.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: difference().
 */
TEST_F(CoefficientMatrixTest, Difference_R) {
  CoefficientMatrix<RESTRICTED> matrix(tr);
  auto diff = matrix.difference();
  EXPECT_EQ(diff(0, 0), 0.0);
  EXPECT_EQ(diff(0, 1), 0.0);
  EXPECT_EQ(diff(1, 0), 0.0);
  EXPECT_EQ(diff(1, 1), 0.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Move Construction
 */
TEST_F(CoefficientMatrixTest, MoveConstruct_R) {
  EXPECT_NE(tr.data(), nullptr);
  CoefficientMatrix<RESTRICTED> result(std::move(tr));
  EXPECT_EQ(result(0, 0), 0.0);
  EXPECT_EQ(result(0, 1), 0.1);
  EXPECT_EQ(result(1, 1), 1.1);
  EXPECT_NE(tr.data(), result.data());
  EXPECT_EQ(result.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Move assignement operator
 */
TEST_F(CoefficientMatrixTest, MoveAssignement_R) {
  CoefficientMatrix<RESTRICTED> result(tBasisController);
  result = std::move(tr);
  EXPECT_EQ(result(0, 0), 0.0);
  EXPECT_EQ(result(0, 1), 0.1);
  EXPECT_EQ(result(1, 1), 1.1);
  EXPECT_NE(tr.data(), result.data());
  EXPECT_EQ(0.0, tr(0, 0));
  EXPECT_EQ(result.getBasisController(), tBasisController);
}

/**
 * @test
 * @brief Tests CoefficientMatrix.h: to cout
 */
TEST_F(CoefficientMatrixTest, Stream_R) {
  std::cout << tr << std::endl;
}

/* ==================
 *    UNRESTRICTED
 * ==================
 */

/**
 * @test
 * @brief Tests CoefficientMatrix.h: Normal construction
 */
TEST_F(CoefficientMatrixTest, BasisConstruct_U) {
  CoefficientMatrix<UNRESTRICTED> matrix(tBasisController);
  // Check for zero initialization
  EXPECT_EQ(matrix.alpha(0, 0), 0.0);
  EXPECT_EQ(matrix.alpha(0, 1), 0.0);
  EXPECT_EQ(matrix.alpha(1, 0), 0.0);
  EXPECT_EQ(matrix.alpha(1, 1), 0.0);
  EXPECT_EQ(matrix.beta(0, 0), 0.0);
  EXPECT_EQ(matrix.beta(1, 0), 0.0);
  EXPECT_EQ(matrix.beta(0, 1), 0.0);
  EXPECT_EQ(matrix.beta(1, 1), 0.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Copy construction.
 */
TEST_F(CoefficientMatrixTest, CopyConstruct_U) {
  CoefficientMatrix<UNRESTRICTED> matrix(tu);
  // Check for zero initialization
  EXPECT_EQ(matrix.alpha(0, 0), 0.0);
  EXPECT_EQ(matrix.alpha(0, 1), 0.1);
  EXPECT_EQ(matrix.alpha(1, 0), 1.0);
  EXPECT_EQ(matrix.alpha(1, 1), 1.1);
  EXPECT_EQ(matrix.beta(0, 0), 0.0);
  EXPECT_EQ(matrix.beta(1, 0), 1.0);
  EXPECT_EQ(matrix.beta(0, 1), 0.1);
  EXPECT_EQ(matrix.beta(1, 1), 1.1);
  EXPECT_EQ(tu.alpha(0, 0), 0.0);
  EXPECT_EQ(tu.alpha(0, 1), 0.1);
  EXPECT_EQ(tu.alpha(1, 0), 1.0);
  EXPECT_EQ(tu.alpha(1, 1), 1.1);
  EXPECT_EQ(tu.beta(0, 0), 0.0);
  EXPECT_EQ(tu.beta(1, 0), 1.0);
  EXPECT_EQ(tu.beta(0, 1), 0.1);
  EXPECT_EQ(tu.beta(1, 1), 1.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: getBasisController(), getNBasisFunctions()
 */
TEST_F(CoefficientMatrixTest, GetBasis_U) {
  CoefficientMatrix<UNRESTRICTED> matrix(tBasisController);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
  EXPECT_EQ(matrix.getNBasisFunctions(), (unsigned int)2);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: total()
 */
TEST_F(CoefficientMatrixTest, Total_U) {
  CoefficientMatrix<UNRESTRICTED> matrix(tu);
  auto total = matrix.total();
  EXPECT_EQ(total(0, 0), 0.0);
  EXPECT_EQ(total(0, 1), 0.2);
  EXPECT_EQ(total(1, 0), 2.0);
  EXPECT_EQ(total(1, 1), 2.2);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: difference().
 */
TEST_F(CoefficientMatrixTest, Difference_U) {
  CoefficientMatrix<UNRESTRICTED> matrix(tu);
  matrix.beta.setZero();
  auto diff = matrix.difference();
  EXPECT_EQ(diff(0, 0), 0.0);
  EXPECT_EQ(diff(0, 1), 0.1);
  EXPECT_EQ(diff(1, 0), 1.0);
  EXPECT_EQ(diff(1, 1), 1.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Move assignement operator
 */
TEST_F(CoefficientMatrixTest, MoveAssignement_U) {
  CoefficientMatrix<UNRESTRICTED> result(tBasisController);
  result = std::move(tu);
  EXPECT_EQ(result.alpha(0, 0), 0.0);
  EXPECT_EQ(result.alpha(0, 1), 0.1);
  EXPECT_EQ(result.alpha(1, 1), 1.1);
  EXPECT_NE(tu.alpha.data(), result.alpha.data());
  EXPECT_EQ(0.0, tu.alpha(0, 0));
  EXPECT_EQ(result.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests CoefficientMatrix.h: Move Construction
 */
TEST_F(CoefficientMatrixTest, MoveConstruct_U) {
  EXPECT_NE(tu.alpha.data(), nullptr);
  CoefficientMatrix<UNRESTRICTED> result(std::move(tu));
  EXPECT_EQ(result.alpha(0, 0), 0.0);
  EXPECT_EQ(result.alpha(0, 1), 0.1);
  EXPECT_EQ(result.alpha(1, 1), 1.1);
  EXPECT_NE(tu.alpha.data(), result.alpha.data());
  EXPECT_EQ(result.getBasisController(), tBasisController);
}

/**
 * @test
 * @brief Tests CoefficientMatrix.h: to cout
 */
TEST_F(CoefficientMatrixTest, Stream_U) {
  std::cout << tu << std::endl;
}
} /* namespace Serenity */
