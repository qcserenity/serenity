/**
 * @file   ScalarProductMatrix_test.cpp
 *
 * @date   Oct 2, 2014
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
#include "math/ScalarProductMatrix.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {
/**
 * @test
 * @brief Tests ScalarProductMatrix.h: creation, filling, deletion
 */
TEST(ScalarProductMatrixTest, CreateFillErase) {
  ScalarProductMatrix testSPM(9, 1);
  Matrix<double> dataOne(9, 1);
  Matrix<double> dataTwo(9, 1);
  Matrix<double> dataThree(9, 1);
  dataOne.fill(1.0);
  dataTwo.fill(2.0);
  dataThree.fill(3.0);
  testSPM.putNewData(dataOne);
  testSPM.putNewData(dataTwo);
  testSPM.putNewData(dataThree);
  EXPECT_EQ(testSPM.getRawData().size(), (unsigned int)3);
  testSPM.eraseAllData();
  EXPECT_EQ(testSPM.getRawData().size(), (unsigned int)0);
}

/**
 * @test
 * @brief Tests ScalarProductMatrix.h: getScalarProductMatrix
 */
TEST(ScalarProductMatrixTest, GetScalarProductMatrix) {
  ScalarProductMatrix testSPM(9, 1);
  Matrix<double> dataOne(9, 1);
  Matrix<double> dataTwo(9, 1);
  Matrix<double> dataThree(9, 1);
  dataOne.fill(1.0);
  dataTwo.fill(2.0);
  dataThree.fill(3.0);
  testSPM.putNewData(dataOne);
  testSPM.putNewData(dataTwo);
  testSPM.putNewData(dataThree);
  auto result = testSPM.getScalarProductMatrix();
  EXPECT_EQ(result(0, 0), 9);
  EXPECT_EQ(result(1, 0), 18);
  EXPECT_EQ(result(0, 1), 18);
  EXPECT_EQ(result(1, 1), 36);
  EXPECT_EQ(result(2, 0), 27);
  EXPECT_EQ(result(0, 2), 27);
  EXPECT_EQ(result(2, 1), 54);
  EXPECT_EQ(result(1, 2), 54);
  EXPECT_EQ(result(2, 2), 81);
}

} /* namespace Serenity */
