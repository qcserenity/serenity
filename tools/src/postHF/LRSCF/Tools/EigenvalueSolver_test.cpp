/**
 * @file   EigenvalueSolver_test.cpp
 *
 * @date   May 30, 2019
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class EigenvalueSolverTest : public ::testing::Test {
 protected:
  EigenvalueSolverTest() {
  }
  virtual ~EigenvalueSolverTest() = default;
};

TEST_F(EigenvalueSolverTest, TDA) {
  // Construct test matrices A+B and A-B
  unsigned nDimension = 250;
  Eigen::MatrixXd A(nDimension, nDimension);
  A.setRandom();
  // Make symmetric
  A = A.transpose() * A;
  // Make diagonal dominant and positive definite
  Eigen::VectorXd diagonalElements = A.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  A *= 1.0e-3;
  A.diagonal() = diagonalElements;

  // Sigma vector calculator
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(1));
    (*sigma)[0] = A * guessVectors[0];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 1000, nEigen, Options::LR_METHOD::TDA,
                          sigmaVectorCalculator);
  Eigen::VectorXd itEigenvalues = solver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(A);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues();

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

TEST_F(EigenvalueSolverTest, RPA) {
  // Construct test matrices A+B and A-B
  unsigned int nDimension = 250;
  Eigen::MatrixXd APB(nDimension, nDimension);
  APB.setRandom();
  // make symmetric
  APB = APB.transpose() * APB;
  // make diagonal dominant and A-B positive definite
  Eigen::VectorXd diagonalElements = APB.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  APB *= 1.0e-3;
  APB.diagonal() = diagonalElements;
  // Assume A-B is diagonal (which is the case when no exact exchange is used)
  Eigen::MatrixXd AMB = diagonalElements.asDiagonal();

  // Sigma vector calculators
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(2));
    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 1000, nEigen, Options::LR_METHOD::TDDFT,
                          sigmaVectorCalculator);
  Eigen::VectorXd itEigenvalues = solver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(2 * nDimension, 2 * nDimension);
  matrix.bottomLeftCorner(nDimension, nDimension) = APB;
  matrix.topRightCorner(nDimension, nDimension) = AMB;
  Eigen::EigenSolver<Eigen::MatrixXd> eigenvalueSolver(matrix, true);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues().real();

  // Sort eigenvalues in ascending order
  std::vector<std::pair<double, unsigned int>> evPair;
  for (unsigned int i = 0; i < 2 * nDimension; ++i) {
    if (exEigenvalues(i) > 0)
      evPair.push_back(std::make_pair(exEigenvalues(i), i));
  }
  std::stable_sort(evPair.begin(), evPair.end());

  // Reassemble eigenvalues and expansion coefficients
  for (unsigned int i = 0; i < nEigen; ++i) {
    exEigenvalues(i) = evPair[i].first;
  }

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

TEST_F(EigenvalueSolverTest, TDASubspaceCollapse) {
  // Construct test matrices A+B and A-B
  unsigned nDimension = 250;
  Eigen::MatrixXd A(nDimension, nDimension);
  A.setRandom();
  // Make symmetric
  A = A.transpose() * A;
  // Make diagonal dominant and positive definite
  Eigen::VectorXd diagonalElements = A.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  A *= 1.0e-3;
  A.diagonal() = diagonalElements;

  // Sigma vector calculator
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(1));
    (*sigma)[0] = A * guessVectors[0];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 60, nEigen, Options::LR_METHOD::TDA,
                          sigmaVectorCalculator);
  Eigen::VectorXd itEigenvalues = solver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(A);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues();

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

TEST_F(EigenvalueSolverTest, TDAInitialGuess) {
  // Construct test matrices A+B and A-B
  unsigned nDimension = 250;
  Eigen::MatrixXd A(nDimension, nDimension);
  A.setRandom();
  // Make symmetric
  A = A.transpose() * A;
  // Make diagonal dominant and positive definite
  Eigen::VectorXd diagonalElements = A.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  A *= 1.0e-3;
  A.diagonal() = diagonalElements;

  // Sigma vector calculator
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(1));
    (*sigma)[0] = A * guessVectors[0];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-3, 100, 1000, nEigen, Options::LR_METHOD::TDA,
                          sigmaVectorCalculator);

  auto guessPtr = std::make_shared<std::vector<Eigen::MatrixXd>>(solver.getEigenvectors());

  EigenvalueSolver guessSolver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 1000, nEigen, Options::LR_METHOD::TDA,
                               sigmaVectorCalculator, guessPtr);
  Eigen::VectorXd itEigenvalues = guessSolver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(A);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues();

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

TEST_F(EigenvalueSolverTest, RPASubspaceCollapse) {
  // Construct test matrices A+B and A-B
  unsigned int nDimension = 250;
  Eigen::MatrixXd APB(nDimension, nDimension);
  APB.setRandom();
  // make symmetric
  APB = APB.transpose() * APB;
  // make diagonal dominant and A-B positive definite
  Eigen::VectorXd diagonalElements = APB.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  APB *= 1.0e-3;
  APB.diagonal() = diagonalElements;
  // Assume A-B is diagonal (which is the case when no exact exchange is used)
  Eigen::MatrixXd AMB = diagonalElements.asDiagonal();

  // Sigma vector calculators
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(2));
    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 60, nEigen, Options::LR_METHOD::TDDFT,
                          sigmaVectorCalculator);
  Eigen::VectorXd itEigenvalues = solver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(2 * nDimension, 2 * nDimension);
  matrix.bottomLeftCorner(nDimension, nDimension) = APB;
  matrix.topRightCorner(nDimension, nDimension) = AMB;
  Eigen::EigenSolver<Eigen::MatrixXd> eigenvalueSolver(matrix, true);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues().real();

  // Sort eigenvalues in ascending order
  std::vector<std::pair<double, unsigned int>> evPair;
  for (unsigned int i = 0; i < 2 * nDimension; ++i) {
    if (exEigenvalues(i) > 0)
      evPair.push_back(std::make_pair(exEigenvalues(i), i));
  }
  std::stable_sort(evPair.begin(), evPair.end());

  // Reassemble eigenvalues and expansion coefficients
  for (unsigned int i = 0; i < nEigen; ++i) {
    exEigenvalues(i) = evPair[i].first;
  }

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

TEST_F(EigenvalueSolverTest, RPAInitialGuess) {
  // Construct test matrices A+B and A-B
  unsigned int nDimension = 250;
  Eigen::MatrixXd APB(nDimension, nDimension);
  APB.setRandom();
  // make symmetric
  APB = APB.transpose() * APB;
  // make diagonal dominant and A-B positive definite
  Eigen::VectorXd diagonalElements = APB.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  APB *= 1.0e-3;
  APB.diagonal() = diagonalElements;
  // Assume A-B is diagonal (which is the case when no exact exchange is used)
  Eigen::MatrixXd AMB = diagonalElements.asDiagonal();

  // Sigma vector calculators
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(2));
    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    return sigma;
  };

  unsigned nEigen = 20;

  // Solve for some eigenvalues
  EigenvalueSolver solver(nDimension, nEigen, diagonalElements, 1.0e-3, 100, 1000, nEigen, Options::LR_METHOD::TDDFT,
                          sigmaVectorCalculator);

  auto guessPtr = std::make_shared<std::vector<Eigen::MatrixXd>>(solver.getEigenvectors());

  EigenvalueSolver guessSolver(nDimension, nEigen, diagonalElements, 1.0e-8, 100, 1000, nEigen,
                               Options::LR_METHOD::TDDFT, sigmaVectorCalculator, guessPtr);

  Eigen::VectorXd itEigenvalues = guessSolver.getEigenvalues();

  // Solve eigenvalue problem using Eigen
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(2 * nDimension, 2 * nDimension);
  matrix.bottomLeftCorner(nDimension, nDimension) = APB;
  matrix.topRightCorner(nDimension, nDimension) = AMB;
  Eigen::EigenSolver<Eigen::MatrixXd> eigenvalueSolver(matrix, true);
  Eigen::VectorXd exEigenvalues = eigenvalueSolver.eigenvalues().real();

  // Sort eigenvalues in ascending order
  std::vector<std::pair<double, unsigned int>> evPair;
  for (unsigned int i = 0; i < 2 * nDimension; ++i) {
    if (exEigenvalues(i) > 0)
      evPair.push_back(std::make_pair(exEigenvalues(i), i));
  }
  std::stable_sort(evPair.begin(), evPair.end());

  // Reassemble eigenvalues and expansion coefficients
  for (unsigned int i = 0; i < nEigen; ++i) {
    exEigenvalues(i) = evPair[i].first;
  }

  // Compare
  for (unsigned int i = 0; i < nEigen; ++i) {
    EXPECT_NEAR(exEigenvalues(i), itEigenvalues(i), 1.0e-7);
  }
}

} /* namespace Serenity */
