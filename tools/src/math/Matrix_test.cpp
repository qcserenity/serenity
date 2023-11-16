/**
 * @file Matrix_test.cpp
 *
 * @date Feb 19, 2014
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
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class MatrixTest
 * @brief Sets everything up for the tests of Matrix.h .
 */
class MatrixTest : public ::testing::Test {
 protected:
  MatrixTest() : tMatrix(4, 3), tMatrix2(3, 2) {
    tMatrix(0, 0) = 1.0;
    tMatrix(0, 1) = 1.0;
    tMatrix(0, 2) = 1.0;
    tMatrix(1, 0) = 2.0;
    tMatrix(1, 1) = 2.0;
    tMatrix(1, 2) = 2.0;
    tMatrix(2, 0) = 3.0;
    tMatrix(2, 1) = 3.0;
    tMatrix(2, 2) = 3.0;
    tMatrix(3, 0) = 4.0;
    tMatrix(3, 1) = 4.0;
    tMatrix(3, 2) = 4.0;
    tMatrix2(0, 0) = 1.0;
    tMatrix2(1, 0) = 1.0;
    tMatrix2(2, 0) = 1.0;
    tMatrix2(0, 1) = 2.0;
    tMatrix2(1, 1) = 2.0;
    tMatrix2(2, 1) = 2.0;
  }

  virtual ~MatrixTest() = default;
  /// The Matrix under test
  Matrix<double> tMatrix;
  Matrix<double> tMatrix2;
};

/**
 * @test
 * @brief Tests Matrix.h: set matrix to matrix (operator =)
 */
TEST_F(MatrixTest, SetMatrix) {
  Matrix<double> result = tMatrix;
  EXPECT_EQ(result(0, 0), 1.0);
  EXPECT_EQ(result(1, 1), 2.0);
  EXPECT_EQ(result(2, 2), 3.0);
  EXPECT_EQ(result(0, 2), 1.0);
  EXPECT_NE(result(2, 2), 10.0);
}

/**
 * @test
 * @brief Tests Matrix.h: multiply matrix with matrix (operator *)
 */
TEST_F(MatrixTest, MultiplyMatrix) {
  Matrix<double> result = (tMatrix * tMatrix2);
  EXPECT_EQ(result(0, 0), 3.0);
  EXPECT_EQ(result(1, 1), 12.0);
  EXPECT_EQ(result(2, 1), 18.0);
  EXPECT_EQ(result(2, 0), 9.0);
  EXPECT_EQ(result(3, 0), 12.0);
  //  EXPECT_EQ(result(0,4), 15.0);
  //  EXPECT_EQ(result(3,4), 60.0);
  EXPECT_NE(result(1, 1), 10.0);
}

TEST_F(MatrixTest, MultiplyEqMatrix) {
  Matrix<double> result = tMatrix;
  result *= tMatrix2;
  EXPECT_EQ(result(0, 0), 3.0);
  EXPECT_EQ(result(1, 1), 12.0);
  EXPECT_EQ(result(2, 1), 18.0);
  EXPECT_EQ(result(2, 0), 9.0);
  EXPECT_EQ(result(3, 0), 12.0);
  //  EXPECT_EQ(result(0,4), 15.0);
  //  EXPECT_EQ(result(3,4), 60.0);
  EXPECT_NE(result(1, 1), 10.0);
}

/**
 * @test
 * @brief Tests Matrix.h: get column
 */
TEST_F(MatrixTest, GetColumn) {
  Matrix<double> result = tMatrix;
  EXPECT_EQ(result.col(1)[0], 1.0);
  EXPECT_EQ(result.col(1)[1], 2.0);
  EXPECT_EQ(result.col(1)[2], 3.0);
  EXPECT_EQ(result.col(1)[3], 4.0);
  EXPECT_NE(result.col(2)[1], 5.0);
}

/**
 * @test
 * @brief Tests Matrix.h: get row
 */
TEST_F(MatrixTest, GetRow) {
  Matrix<double> result = tMatrix;
  EXPECT_EQ(result.row(2)[0], 3.0);
  EXPECT_EQ(result.row(2)[1], 3.0);
  EXPECT_EQ(result.row(2)[2], 3.0);
  EXPECT_NE(result.row(2)[1], 5.0);
}

/**
 * @test
 * @brief Tests Matrix.h: multiply matrix with scalar (operator *)
 */
TEST_F(MatrixTest, MultiplyScalar) {
  Matrix<double> result = tMatrix * 5.0;
  EXPECT_EQ(result(0, 0), 5.0);
  EXPECT_EQ(result(1, 1), 10.0);
  EXPECT_EQ(result(2, 2), 15.0);
  EXPECT_EQ(result(0, 2), 5.0);
  EXPECT_NE(result(2, 2), 156.0);

  result = 5.0 * tMatrix;
  EXPECT_EQ(result(0, 0), 5.0);
  EXPECT_EQ(result(1, 1), 10.0);
  EXPECT_EQ(result(2, 2), 15.0);
  EXPECT_EQ(result(0, 2), 5.0);
  EXPECT_EQ(result(3, 2), 20.0);
  EXPECT_NE(result(2, 2), 156.0);
}

/**
 * @test
 * @brief Tests Matrix.h: multiply matrix with scalar (operator *=)
 */
TEST_F(MatrixTest, MultiplyEqScalar) {
  Matrix<double> result = tMatrix;
  result *= 5.0;
  EXPECT_EQ(result(0, 0), 5.0);
  EXPECT_EQ(result(1, 1), 10.0);
  EXPECT_EQ(result(2, 2), 15.0);
  EXPECT_EQ(result(0, 2), 5.0);
  EXPECT_EQ(result(3, 2), 20.0);
  EXPECT_NE(result(2, 2), 156.0);
}

/**
 * @test
 * @brief Tests Matrix.h: addition of matrix and matrix (operator +)
 */
TEST_F(MatrixTest, AddMatrix) {
  Matrix<double> result = tMatrix + tMatrix;
  EXPECT_EQ(result(0, 0), 2.0);
  EXPECT_EQ(result(1, 1), 4.0);
  EXPECT_EQ(result(2, 2), 6.0);
  EXPECT_EQ(result(0, 2), 2.0);
  EXPECT_EQ(result(3, 2), 8.0);
  EXPECT_NE(result(2, 2), 10.0);
}

/**
 * @test
 * @brief Tests Matrix.h: addition of matrix and matrix (operator +=)
 */
TEST_F(MatrixTest, AddEqMatrix) {
  Matrix<double> result = tMatrix;
  result += tMatrix;
  EXPECT_EQ(result(0, 0), 2.0);
  EXPECT_EQ(result(1, 1), 4.0);
  EXPECT_EQ(result(2, 2), 6.0);
  EXPECT_EQ(result(0, 2), 2.0);
  EXPECT_EQ(result(3, 2), 8.0);
  EXPECT_NE(result(2, 2), 10.0);
}

/**
 * @test
 * @brief Tests Matrix.h: subtraction of matrix and matrix (operator -)
 */
TEST_F(MatrixTest, SubtractMatrix) {
  Matrix<double> result = tMatrix - tMatrix;
  EXPECT_EQ(result(0, 0), 0.0);
  EXPECT_EQ(result(1, 1), 0.0);
  EXPECT_EQ(result(2, 2), 0.0);
  EXPECT_EQ(result(0, 2), 0.0);
  EXPECT_EQ(result(3, 2), 0.0);
  EXPECT_NE(result(2, 2), 10.0);
}

/**
 * @test
 * @brief Tests Matrix.h: subtraction of matrix and matrix (operator -=)
 */
TEST_F(MatrixTest, SubtractEqMatrix) {
  Matrix<double> result = tMatrix;
  result -= tMatrix;
  EXPECT_EQ(result(0, 0), 0.0);
  EXPECT_EQ(result(1, 1), 0.0);
  EXPECT_EQ(result(2, 2), 0.0);
  EXPECT_EQ(result(0, 2), 0.0);
  EXPECT_EQ(result(3, 2), 0.0);
  EXPECT_NE(result(2, 2), 10.0);
}

} // namespace Serenity
