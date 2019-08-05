/**
 * @file IterativeEigenvalueSolver_test.cpp
 *
 * @date Oct 17, 2016
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
#include "postHF/LRSCF/Solver/IterativeEigenvalueSolver.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>



namespace Serenity {

class IterativeEigenvalueSolverTest: public :: testing::Test {
protected:
  IterativeEigenvalueSolverTest() {}
  virtual ~IterativeEigenvalueSolverTest() = default;
};




/**
 * @brief Iterative eigenvalue solver class for tests. This is necessary, since the objects tested here are
 *        protected and can only be accessed from inheriting classes.
 */
class TestSolver: public IterativeEigenvalueSolver {

public:
  TestSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations);

  virtual ~TestSolver() = default;

  std::vector<Eigen::MatrixXd >& getGuessVectors();

  void computeTrials(
      Eigen::MatrixXd& trialVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::MatrixXd& guessVectors);

  void computeResidual(
      Eigen::MatrixXd& residualVectors,
      Eigen::MatrixXd& sigmaVectors,
      Eigen::MatrixXd& guessVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::VectorXd& subspaceEigenvalues);

protected:
  void iterate() override final;
  void postProcessing() override final;
};

  TestSolver::TestSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations):
        IterativeEigenvalueSolver(nDimension,nEigen,diagonalElements,convergenceCriterion,nMaxIterations) {
    _guessVectors.resize(1);
  }

  std::vector<Eigen::MatrixXd >&TestSolver::getGuessVectors() {
    computeInitialGuessVectors(_guessVectors[0],_diagonalElements);
    return _guessVectors;
  }

  void TestSolver::computeTrials(
      Eigen::MatrixXd& trialVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::MatrixXd& guessVectors) {
    computeTrialVectors(trialVectors,expansionCoefficients,guessVectors);
  }

  void TestSolver::computeResidual(
      Eigen::MatrixXd& residualVectors,
      Eigen::MatrixXd& sigmaVectors,
      Eigen::MatrixXd& guessVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::VectorXd& subspaceEigenvalues) {
    computeResidualVectors(residualVectors,sigmaVectors,guessVectors,expansionCoefficients,subspaceEigenvalues);
  }

  void TestSolver::iterate() {return;}

  void TestSolver::postProcessing(){return;}






TEST_F(IterativeEigenvalueSolverTest, IterativeEigenvalueSolver) {
  /*
   * The solve function is already tested by actual implementations of iterative eigenvalue solvers.
   * What is tested here, are the helper functions which come along with the iterative eigenvalue solver
   * class.
   */

  //set up primitive diagonal element vector
  Eigen::VectorXd diagonal(5);
  diagonal(0) = 0.3;
  diagonal(1) = 0.2;
  diagonal(2) = 0.1;
  diagonal(3) = 0.4;
  diagonal(4) = 0.5;
  //get instance of test eigenvalue solver
  TestSolver testSolver(
      diagonal.rows(),
      diagonal.rows(),
      diagonal,
      1.0e-5,
      100);

  /*
   * Test unit guess vectors
   */

  std::vector<Eigen::MatrixXd >& guessVector = testSolver.getGuessVectors();
  //compare
  EXPECT_EQ(guessVector[0](0),0.0);
  EXPECT_EQ(guessVector[0](1),0.0);
  EXPECT_EQ(guessVector[0](2),1.0);
  EXPECT_EQ(guessVector[0](3),0.0);
  EXPECT_EQ(guessVector[0](4),0.0);

  /*
   * Test trial vectors: trialVectors = guessVectors * expansionCoefficients
   */
  unsigned int nDimension = diagonal.rows();
  Eigen::MatrixXd randomMatrix1(nDimension,nDimension);
  Eigen::MatrixXd randomMatrix2(nDimension,nDimension);
  randomMatrix1.setRandom();
  randomMatrix2.setRandom();
  Eigen::MatrixXd trials1 = randomMatrix1 * randomMatrix2;
  Eigen::MatrixXd trials2;
  testSolver.computeTrials(trials2,randomMatrix2,randomMatrix1);
  //compare
  Eigen::MatrixXd thisShouldBeZero1 = trials1 - trials2;
  for (unsigned int i = 0; i < diagonal.rows(); ++i) {
    for (unsigned int j = 0; j < diagonal.rows(); ++j) {
      EXPECT_NEAR(thisShouldBeZero1(i,j),0.0,1.0e-7);
    }
  }


  /*
   * Test residual vectors (compare with hard coded example)
   */

  Eigen::MatrixXd expansionCoefficients(nDimension,nDimension);
  expansionCoefficients.setRandom();
  Eigen::MatrixXd sigmaVectors(nDimension,nDimension);
  sigmaVectors.setRandom();
  Eigen::MatrixXd guessVectors(nDimension,nDimension);
  guessVectors.setRandom();
  //calculate residual of first vector
  Eigen::MatrixXd differenceMatrix(nDimension,nDimension);
  differenceMatrix.setZero();
  for (unsigned int i = 0; i < nDimension; ++i) {
    differenceMatrix.col(i) = sigmaVectors.col(i) - diagonal(0) * guessVectors.col(i);
  }
  Eigen::VectorXd residual(nDimension);
  residual.setZero();
  for (unsigned int i = 0; i < nDimension; ++i) {
    residual += expansionCoefficients(i,0) * differenceMatrix.col(i);
  }
  Eigen::MatrixXd residuals;
  testSolver.computeResidual(residuals,sigmaVectors,guessVectors,expansionCoefficients,diagonal);
  //compare
  Eigen::VectorXd thisShouldBeZero2 = residual - residuals.col(0);
  for (unsigned int i = 0; i < nDimension; ++i) {
    EXPECT_NEAR(thisShouldBeZero2(i),0.0,1.0e-7);
  }

}

} /* namespace Serenity */
