/**
 * @file DavidsonSolver_test.cpp
 *
 * @date Aug 31, 2016
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
#include "postHF/LRSCF/Solver/DavidsonSolver.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {
/**
 * @class DavidsonSolverTest
 * @brief Test Davidson eigenvalue solver by calculating eigenpairs of a symmetric diagonal
 *        dominant random matrix and comparing with the exact result.
 */
class DavidsonSolverTest : public ::testing::Test {
protected:
  DavidsonSolverTest() {}
  virtual ~DavidsonSolverTest() = default;
};

TEST_F(DavidsonSolverTest,DavidsonEigenvalueSolver) {
  //set up random matrix for test diagonalization
  unsigned int nDimension = 100;
  Eigen::MatrixXd randomMatrix(nDimension,nDimension);
  randomMatrix.setRandom();
  //Work's only for symmetric matrices
  randomMatrix = randomMatrix.transpose() * randomMatrix;
  Eigen::VectorXd diagonalElements = randomMatrix.diagonal();
  //Work's only for diagonal dominant matrices
  randomMatrix *= 1.0e-3;
  diagonalElements *= 1.0e-1;
  randomMatrix.diagonal() = diagonalElements;
  //Sigma vector calculator for Davidson solver
  auto sigmaVectorCalculator = [&] (
      Eigen::MatrixXd& guessVectors) {
    Eigen::MatrixXd sigmaVectors;
    sigmaVectors = randomMatrix * guessVectors;
    return sigmaVectors;
  };
  //get 5 lowest eigenpairs
  DavidsonSolver davidsonSolver(
      nDimension,
      5,
      diagonalElements,
      1.0e-8,
      999,
      sigmaVectorCalculator,
      50);
  davidsonSolver.solve();
  Eigen::VectorXd eigenvalues = davidsonSolver.getEigenvalues();
  Eigen::MatrixXd eigenvectors = davidsonSolver.getEigenvectors()[0];
  //exactly diagonlaize randomMatrix using eigen's qr algorithm
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(nDimension);
  eigenvalueSolver.compute(randomMatrix);
  Eigen::VectorXd exactEigenvalues = eigenvalueSolver.eigenvalues();
  Eigen::MatrixXd exactEigenvectors = eigenvalueSolver.eigenvectors();
  //Compare eigenvalues
  for (unsigned int i = 0; i < eigenvalues.rows(); ++i) {
    EXPECT_NEAR(eigenvalues(i),exactEigenvalues(i),1.0e-7);
  }
  //Compare eigenvectors
  for (unsigned int i = 0; i < eigenvectors.cols(); ++i) {
    for (unsigned int j = 0; j < nDimension; ++j) {
      //take absolute values, eigenvectors are not unique in sign
      double thisShouldBeZero = fabs(eigenvectors(j,i)) - fabs(exactEigenvectors(j,i));
      EXPECT_NEAR(thisShouldBeZero,0.0,1.0e-7);
    }
  }
}



} /* namespace Serenity */
