/**
 * @file CholeskyDecomposer_test.cpp
 *
 * @date Jan 25, 2017
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
#include "math/linearAlgebra/CholeskyDecomposer.h"
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/Sparse>
#include <gtest/gtest.h>


namespace Serenity {

class CholeskyDecomposerTest: public::testing::Test {
protected:
  CholeskyDecomposerTest(){}
  virtual ~CholeskyDecomposerTest() = default;
};

TEST_F(CholeskyDecomposerTest,CholeskyDecomposition) {
  //build random positive definite matrix
  Eigen::MatrixXd randomMatrix(100,100);
  randomMatrix.setRandom();
  randomMatrix = randomMatrix * randomMatrix.transpose();
  //get diagonal of randomMatrix
  Eigen::VectorXd diagonal = randomMatrix.diagonal();
  //function to calculate matrix cols of randomMatrix
  auto columnCalculator = [&] (
      std::vector<unsigned int>& qualifiedSet) {
    auto colMat = std::unique_ptr<Eigen::MatrixXd>(
        new Eigen::MatrixXd(diagonal.rows(),qualifiedSet.size()));
    for (unsigned int p = 0; p < diagonal.rows(); ++p) {
      for (unsigned int q = 0; q < qualifiedSet.size(); ++q) {
        (*colMat)(p,q) = randomMatrix(p,qualifiedSet[q]);
      }
    }
    return colMat;
  };
  //choose large span factor so that decomposer has to decompose many blocks
  CholeskyDecomposer decomposer(diagonal,columnCalculator,1.0e-8,10,0.1);
  auto choleskyVectors = decomposer.getCholeskyVectors();
  //reconstruct matrix from Cholesky vectors
  Eigen::SparseMatrix<double>  testMatrix = (*choleskyVectors) * (*choleskyVectors).transpose();
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

} /* namespace Serenity */
