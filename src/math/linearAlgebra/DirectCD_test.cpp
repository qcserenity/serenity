/**
 * @file DirectCD_test.cpp
 *
 * @date Jun 10, 2016
 * @author Michael Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Serenity Internal Headers */
#include "math/linearAlgebra/DirectCD.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/Sparse>
#include <gtest/gtest.h>


namespace Serenity {
/**
 * @class DirectCDTest
 * @brief Tests direct pivoted Cholesky decomposition by decomposing a random
 *        positive definite matrix and comparing it with it's reconstructed
 *        matrix.
 */
class DirectCDTest : public ::testing::Test {
protected:
  DirectCDTest(){}
  virtual ~DirectCDTest() = default;
};

TEST_F(DirectCDTest,CholeskyDecomposition) {
  //build random positive definite matrix
  Matrix<double> randomMatrix(20,20);
  randomMatrix.setRandom();
  randomMatrix = randomMatrix * randomMatrix.transpose();
  //get diagonal of randomMatrix
  VectorXd diagonal = randomMatrix.diagonal();
  //function to calculate matrix cols of randomMatrix
  auto columnCalculator = [&] (VectorXi& qualifiedSet) {
    Matrix<double> colMat(randomMatrix.rows(),qualifiedSet.rows());
    for (unsigned int i=0; i < qualifiedSet.rows(); ++i) {
      colMat.col(i)=randomMatrix.col(qualifiedSet(i));
    }
    return colMat;
  };
  //choose large span factor so that directCD has to decompose many blocks
  DirectCD directCD(diagonal,columnCalculator,1.0e-8,0.9);
  directCD.decompose();
  SparseMatrix<double> choleskyVectors = directCD.getCholeskyVectors();
  //reconstruct matrix from Cholesky vectors
  SparseMatrix<double> testMatrix = choleskyVectors * choleskyVectors.transpose();
  //Compare. It is not necessary to compare every coefficient, one could e.g. also
  //calculate the Frobenius norm of the difference between randomMatrix and testMatrix, which
  //should be close to zero...
  for (unsigned int i = 0; i < randomMatrix.rows();++i) {
    for (unsigned int j = 0; j < randomMatrix.cols();++j){
      double val = randomMatrix(i,j);
      double testVal = testMatrix.coeffRef(i,j);
      EXPECT_NEAR(testVal,val,1.0e-7);
    }
  }
}



}


