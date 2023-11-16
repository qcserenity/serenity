/**
 * @file   DensityMatrix_test.cpp
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
#include "data/matrices/DensityMatrix.h"
#include "testsupply/BasisController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class DensityMatrixTest : public ::testing::Test {
 protected:
  DensityMatrixTest()
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
  virtual ~DensityMatrixTest() = default;
  // The basis in which it is defined
  std::shared_ptr<BasisController> tBasisController;
  // The object under test
  DensityMatrix<RESTRICTED> tr;
  DensityMatrix<UNRESTRICTED> tu;
};
/**
 * @test
 * @brief Tests DensityMatrix.h: Normal construction.
 */
TEST_F(DensityMatrixTest, BasisConstruct_R) {
  DensityMatrix<RESTRICTED> matrix(tBasisController);
  // Check for zero initialization
  EXPECT_EQ(matrix(0, 0), 0.0);
  EXPECT_EQ(matrix(0, 1), 0.0);
  EXPECT_EQ(matrix(1, 0), 0.0);
  EXPECT_EQ(matrix(1, 1), 0.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Copy construction.
 */
TEST_F(DensityMatrixTest, CopyConstruct_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator= MatrixXd.
 */
TEST_F(DensityMatrixTest, Eq_MatrixXd_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator= MiB.
 */
TEST_F(DensityMatrixTest, Eq_MatInB_3_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator= DensityMatrix.
 */
TEST_F(DensityMatrixTest, Eq_MatInB_1_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  DensityMatrix<RESTRICTED> test(tBasisController);
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
 * @brief Tests DensityMatrix.h: operator= DensityMatrix.
 */
TEST_F(DensityMatrixTest, Eq_Ref_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator= DensityMatrix.
 */
TEST_F(DensityMatrixTest, Eq_MatInB_2_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  auto bas2 = BasisController__TEST_SUPPLY::getBasisController(TEST_BASIS_CONTROLLERS::SMALL_MIXED);
  DensityMatrix<RESTRICTED> test(bas2);
  EXPECT_THROW(matrix.operator=(test), SerenityError);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: operator= Eigen3 computed expression.
 */
TEST_F(DensityMatrixTest, Eq_MatProd_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator+= MatrixXd.
 */
TEST_F(DensityMatrixTest, PE_MatrixXd_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator-= MatrixXd.
 */
TEST_F(DensityMatrixTest, ME_MatrixXd_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
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
 * @brief Tests DensityMatrix.h: operator+= DensityMatrix.
 */
TEST_F(DensityMatrixTest, PE_MatInB_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  DensityMatrix<RESTRICTED> test(tBasisController);
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
 * @brief Tests DensityMatrix.h: operator-= DensityMatrix.
 */
TEST_F(DensityMatrixTest, ME_MatInB_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  DensityMatrix<RESTRICTED> test(tBasisController);
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
 * @brief Tests DensityMatrix.h: operator* DensityMatrix.
 */
TEST_F(DensityMatrixTest, Multiply_MatInB_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  DensityMatrix<RESTRICTED> test(tBasisController);
  test(0, 0) = 1.0;
  test(1, 0) = 2.0;
  test(0, 1) = 2.0;
  test(1, 1) = 3.0;
  auto result = matrix * test;
  EXPECT_EQ(result.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: operator* mixed.
 */
TEST_F(DensityMatrixTest, Multiply_Mixed_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  DensityMatrix<RESTRICTED> test1(tBasisController);
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
 * @brief Tests DensityMatrix.h: getBasisController(), getNBasisFunctions()
 */
TEST_F(DensityMatrixTest, GetBasis_R) {
  DensityMatrix<RESTRICTED> matrix(tBasisController);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
  EXPECT_EQ(matrix.getNBasisFunctions(), (unsigned int)2);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: total()
 */
TEST_F(DensityMatrixTest, Total_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  auto total = matrix.total();
  EXPECT_EQ(total(0, 0), 0.0);
  EXPECT_EQ(total(0, 1), 0.1);
  EXPECT_EQ(total(1, 0), 1.0);
  EXPECT_EQ(total(1, 1), 1.1);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: difference().
 */
TEST_F(DensityMatrixTest, Difference_R) {
  DensityMatrix<RESTRICTED> matrix(tr);
  auto diff = matrix.difference();
  EXPECT_EQ(diff(0, 0), 0.0);
  EXPECT_EQ(diff(0, 1), 0.0);
  EXPECT_EQ(diff(1, 0), 0.0);
  EXPECT_EQ(diff(1, 1), 0.0);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Move Construction
 */
TEST_F(DensityMatrixTest, MoveConstruct_R) {
  EXPECT_NE(tr.data(), nullptr);
  DensityMatrix<RESTRICTED> result(std::move(tr));
  EXPECT_EQ(result(0, 0), 0.0);
  EXPECT_EQ(result(0, 1), 0.1);
  EXPECT_EQ(result(1, 1), 1.1);
  EXPECT_NE(tr.data(), result.data());
  EXPECT_EQ(result.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: Move assignement operator
 */
TEST_F(DensityMatrixTest, MoveAssignement_R) {
  DensityMatrix<RESTRICTED> result(tBasisController);
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
 * @brief Tests DensityMatrix.h: to cout
 */
TEST_F(DensityMatrixTest, Stream_R) {
  std::cout << tr << std::endl;
}

/* ==================
 *    UNRESTRICTED
 * ==================
 */

/**
 * @test
 * @brief Tests DensityMatrix.h: Normal construction
 */
TEST_F(DensityMatrixTest, BasisConstruct_U) {
  DensityMatrix<UNRESTRICTED> matrix(tBasisController);
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
 * @brief Tests DensityMatrix.h: Copy construction.
 */
TEST_F(DensityMatrixTest, CopyConstruct_U) {
  DensityMatrix<UNRESTRICTED> matrix(tu);
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
 * @brief Tests DensityMatrix.h: getBasisController(), getNBasisFunctions()
 */
TEST_F(DensityMatrixTest, GetBasis_U) {
  DensityMatrix<UNRESTRICTED> matrix(tBasisController);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
  EXPECT_EQ(matrix.getNBasisFunctions(), (unsigned int)2);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: total()
 */
TEST_F(DensityMatrixTest, Total_U) {
  DensityMatrix<UNRESTRICTED> matrix(tu);
  auto total = matrix.total();
  EXPECT_EQ(total(0, 0), 0.0);
  EXPECT_EQ(total(0, 1), 0.2);
  EXPECT_EQ(total(1, 0), 2.0);
  EXPECT_EQ(total(1, 1), 2.2);
  EXPECT_EQ(matrix.getBasisController(), tBasisController);
}
/**
 * @test
 * @brief Tests DensityMatrix.h: difference().
 */
TEST_F(DensityMatrixTest, Difference_U) {
  DensityMatrix<UNRESTRICTED> matrix(tu);
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
 * @brief Tests DensityMatrix.h: Move assignement operator
 */
TEST_F(DensityMatrixTest, MoveAssignement_U) {
  DensityMatrix<UNRESTRICTED> result(tBasisController);
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
 * @brief Tests DensityMatrix.h: Move Construction
 */
TEST_F(DensityMatrixTest, MoveConstruct_U) {
  EXPECT_NE(tu.alpha.data(), nullptr);
  DensityMatrix<UNRESTRICTED> result(std::move(tu));
  EXPECT_EQ(result.alpha(0, 0), 0.0);
  EXPECT_EQ(result.alpha(0, 1), 0.1);
  EXPECT_EQ(result.alpha(1, 1), 1.1);
  EXPECT_NE(tu.alpha.data(), result.alpha.data());
  EXPECT_EQ(result.getBasisController(), tBasisController);
}

/**
 * @test
 * @brief Tests DensityMatrix.h: to cout
 */
TEST_F(DensityMatrixTest, Stream_U) {
  std::cout << tu << std::endl;
}

} /* namespace Serenity */
