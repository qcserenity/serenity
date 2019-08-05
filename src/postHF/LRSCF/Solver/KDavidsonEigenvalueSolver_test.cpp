/**
 * @file   KDavidsonEigenvalueSolver_test.cpp
 *
 * @date   Mar 4, 2017
 * @author M. Boeckers
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
#include "postHF/LRSCF/Solver/KDavidsonEigenvalueSolver.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {

using namespace Eigen;
/**
 * Test class for KDavidsonEigenvalueSolver
 */
class KDavidsonSolverTest: public ::testing::Test {
protected:
  KDavidsonSolverTest(){}
  virtual ~KDavidsonSolverTest() = default;
};

TEST_F(KDavidsonSolverTest,KDavidsonEigenvalueSolver) {
  //
  //construct test matrices A+B and A-B
  //
  unsigned int nDimension = 100;
  Eigen::MatrixXd APB(nDimension,nDimension);
  APB.setRandom();
  //make symmetric
  APB = APB.transpose() * APB;
  //make diagonal dominant and A-B positive definite
  Eigen::VectorXd diagonalElements = 0.1 * APB.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements= diagonalElements.array().sqrt();
  APB *= 1.0e-3;
  APB.diagonal() = diagonalElements;
  //Assume A-B is diagonal (which is the case when no exact exchange is used)
  Eigen::MatrixXd AMB = diagonalElements.asDiagonal();


  //
  //Sigma vector calculators
  //
  auto sigmaVectorCalculator_APB = [&] (Eigen::MatrixXd  guessVectors) {
    return APB * guessVectors;
  };
  auto sigmaVectorCalculator_AMB = [&] (Eigen::MatrixXd  guessVectors) {
    return AMB * guessVectors;
  };
  //
  //Solve for some eigenvalues
  //
  KDavidsonEigenvalueSolver kDavidsonSolver(
      nDimension,
      10,
      diagonalElements,
      1.0e-7,
      100,
      sigmaVectorCalculator_APB,
      sigmaVectorCalculator_AMB,
      150);
  kDavidsonSolver.solve();
  Eigen::VectorXd kEigenvalues = kDavidsonSolver.getEigenvalues();

  //
  //Solve (symmetric) eigenvalue problem
  //

  //calculate sqrt of matrix A-B to construct symmetric eigenvalue problem
  Eigen::VectorXd temp = diagonalElements.array().sqrt();
  Eigen::MatrixXd sqrtAMB = temp.asDiagonal();

  Eigen::MatrixXd symMatrix = sqrtAMB * APB * sqrtAMB;
  //solve eigenvalue problem using eigens qr algorithm
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(nDimension);
  eigenvalueSolver.compute(symMatrix);
  auto exactEigenvalues = eigenvalueSolver.eigenvalues();

  exactEigenvalues = exactEigenvalues.array().sqrt();

  //
  //Compare
  //

  //compare eigenvalues
  for (unsigned int i = 0; i < kEigenvalues.rows(); ++i) {
    EXPECT_NEAR(exactEigenvalues(i),kEigenvalues(i),1.0e-7);
  }

}


}



