/**
 * @file   ResponseSolver_test.cpp
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
#include "postHF/LRSCF/Tools/ResponseSolver.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

class ResponseSolverTest : public ::testing::Test {
 protected:
  ResponseSolverTest() {
  }
  virtual ~ResponseSolverTest() = default;
};

TEST_F(ResponseSolverTest, Undamped) {
  // construct test matrices A+B and A-B
  unsigned int nDim = 250;
  Eigen::MatrixXd APB(nDim, nDim);
  Eigen::MatrixXd dipE(nDim, 3);
  dipE.setRandom();
  dipE *= 1000;
  std::vector<Eigen::MatrixXd> ppmq(2, 2 * dipE);
  ppmq[1].setZero();

  double damping = 0.;

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
    std::unique_ptr<std::vector<Eigen::MatrixXd>> sigma =
        std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(2));
    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    return sigma;
  };

  std::vector<double> frequencies = {1};

  // Build response solver
  auto responseSolver = ResponseSolver(diagonalElements, 1e-6, 100, 100000, frequencies, damping, ppmq, sigmaVectorCalculator);

  // Solutionvectors
  auto itVectors = responseSolver.getEigenvectors();

  // Setup 'exact' lhs-matrix to feed eigen3 with
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(2 * nDim, 2 * nDim);
  lhs.block(0 * nDim, 0 * nDim, nDim, nDim) = APB;
  lhs.block(1 * nDim, 1 * nDim, nDim, nDim) = AMB;

  lhs.block(1 * nDim, 0 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(0 * nDim, 1 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);

  // Same for rhs-matrix
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2 * nDim, 3);
  rhs.block(0 * nDim, 0, nDim, 3) = 2 * dipE;
  rhs.block(1 * nDim, 0, nDim, 3) = Eigen::MatrixXd::Zero(nDim, 3);

  // Determine solution vectors using eigen
  Eigen::MatrixXd solVectors = lhs.colPivHouseholderQr().solve(rhs);

  // Build difference matrix
  Eigen::MatrixXd diffXpY = solVectors.middleRows(0 * nDim, nDim) - itVectors[0];
  Eigen::MatrixXd diffXmY = solVectors.middleRows(1 * nDim, nDim) - itVectors[1];

  // Compare
  EXPECT_LT(diffXpY.array().cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT(diffXmY.array().cwiseAbs().maxCoeff(), 1e-10);
}

TEST_F(ResponseSolverTest, Undamped_Multifreq) {
  // construct test matrices A+B and A-B
  unsigned int nDim = 250;
  Eigen::MatrixXd APB(nDim, nDim);
  Eigen::MatrixXd dipE(nDim, 3);
  dipE.setRandom();
  dipE *= 1000;
  std::vector<Eigen::MatrixXd> ppmq(2, 2 * dipE);
  ppmq[1].setZero();

  double damping = 0.;

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
    std::unique_ptr<std::vector<Eigen::MatrixXd>> sigma =
        std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(2));
    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    return sigma;
  };

  // Solve for 400 frequencies simultaneously
  std::vector<double> frequencies;
  for (double freq = 0; freq <= 10; freq += 0.025) {
    frequencies.push_back(freq);
  }

  // Build response solver
  auto responseSolver = ResponseSolver(diagonalElements, 1e-6, 100, 100000, frequencies, damping, ppmq, sigmaVectorCalculator);

  // Solutionvectors
  auto itVectors = responseSolver.getEigenvectors();

  // Setup 'exact' lhs-matrix to feed eigen3 with
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(2 * nDim, 2 * nDim);
  lhs.block(0 * nDim, 0 * nDim, nDim, nDim) = APB;
  lhs.block(1 * nDim, 1 * nDim, nDim, nDim) = AMB;

  // Same for rhs-matrix
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2 * nDim, 3);
  rhs.block(0 * nDim, 0, nDim, 3) = 2 * dipE;
  rhs.block(1 * nDim, 0, nDim, 3) = Eigen::MatrixXd::Zero(nDim, 3);

  Eigen::MatrixXd diffXpY, diffXmY, solVectors;
  // Loop over all frequencies
  for (unsigned i = 0; i < frequencies.size(); ++i) {
    lhs.block(1 * nDim, 0 * nDim, nDim, nDim) = -frequencies[i] * Eigen::MatrixXd::Identity(nDim, nDim);
    lhs.block(0 * nDim, 1 * nDim, nDim, nDim) = -frequencies[i] * Eigen::MatrixXd::Identity(nDim, nDim);

    // Determine solution vectors using eigen
    solVectors = lhs.colPivHouseholderQr().solve(rhs);

    diffXpY = solVectors.middleRows(0 * nDim, nDim) - itVectors[0].middleCols(3 * i, 3);
    diffXmY = solVectors.middleRows(1 * nDim, nDim) - itVectors[1].middleCols(3 * i, 3);

    // Compare
    EXPECT_LT(diffXpY.array().cwiseAbs().maxCoeff(), 1e-10);
    EXPECT_LT(diffXmY.array().cwiseAbs().maxCoeff(), 1e-10);
  }
}

TEST_F(ResponseSolverTest, Damped) {
  // construct test matrices A+B and A-B
  unsigned int nDim = 250;
  Eigen::MatrixXd APB(nDim, nDim);
  Eigen::MatrixXd dipE(nDim, 3);
  dipE.setRandom();
  dipE *= 1000;
  std::vector<Eigen::MatrixXd> ppmq(4, 2 * dipE);
  ppmq[1].setZero();
  ppmq[2].setZero();
  ppmq[3].setZero();

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
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(4));

    (*sigma)[0] = APB * guessVectors[0];
    (*sigma)[1] = AMB * guessVectors[1];
    (*sigma)[2] = APB * guessVectors[2];
    (*sigma)[3] = AMB * guessVectors[3];

    return sigma;
  };

  std::vector<double> frequencies = {1.0};
  double damping = 1.0;

  // Build response solver
  auto responseSolver = ResponseSolver(diagonalElements, 1e-6, 100, 100000, frequencies, damping, ppmq, sigmaVectorCalculator);

  // Solutionvectors
  auto itVectors = responseSolver.getEigenvectors();

  // Setup 'exact' lhs-matrix to feed eigen3 with, only for one frequency right now
  // ToDo: Test functionality for multiple frequencies simultaneously
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(4 * nDim, 4 * nDim);
  lhs.block(0 * nDim, 0 * nDim, nDim, nDim) = APB;
  lhs.block(1 * nDim, 1 * nDim, nDim, nDim) = AMB;
  lhs.block(2 * nDim, 2 * nDim, nDim, nDim) = APB;
  lhs.block(3 * nDim, 3 * nDim, nDim, nDim) = AMB;

  lhs.block(0 * nDim, 3 * nDim, nDim, nDim) = damping * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(1 * nDim, 2 * nDim, nDim, nDim) = damping * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(2 * nDim, 1 * nDim, nDim, nDim) = -damping * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(3 * nDim, 0 * nDim, nDim, nDim) = -damping * Eigen::MatrixXd::Identity(nDim, nDim);

  lhs.block(1 * nDim, 0 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(0 * nDim, 1 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(2 * nDim, 3 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);
  lhs.block(3 * nDim, 2 * nDim, nDim, nDim) = -frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);

  // Same for rhs-matrix
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(4 * nDim, 3);
  rhs.block(0 * nDim, 0, nDim, 3) = 2 * dipE;
  rhs.block(1 * nDim, 0, nDim, 3) = Eigen::MatrixXd::Zero(nDim, 3);
  rhs.block(2 * nDim, 0, nDim, 3) = Eigen::MatrixXd::Zero(nDim, 3);
  rhs.block(3 * nDim, 0, nDim, 3) = Eigen::MatrixXd::Zero(nDim, 3);

  // Determine solution vectors using eigen
  Eigen::MatrixXd solVectors = lhs.colPivHouseholderQr().solve(rhs);

  // Build difference matrix
  Eigen::MatrixXd diffXpYr = solVectors.middleRows(0 * nDim, nDim) - itVectors[0];
  Eigen::MatrixXd diffXmYr = solVectors.middleRows(1 * nDim, nDim) - itVectors[1];
  Eigen::MatrixXd diffXpYi = solVectors.middleRows(2 * nDim, nDim) - itVectors[2];
  Eigen::MatrixXd diffXmYi = solVectors.middleRows(3 * nDim, nDim) - itVectors[3];

  // Compare
  EXPECT_LT(diffXpYr.array().cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT(diffXmYr.array().cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT(diffXpYi.array().cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_LT(diffXmYi.array().cwiseAbs().maxCoeff(), 1e-10);
}

TEST_F(ResponseSolverTest, Single_Set) {
  // construct test matrices A+B and A-B
  unsigned int nDim = 250;
  Eigen::MatrixXd APB(nDim, nDim);
  Eigen::MatrixXd dipE(nDim, 3);
  dipE.setRandom();
  dipE *= 1000;
  std::vector<Eigen::MatrixXd> ppmq(1, 2 * dipE);

  double damping = 0.;

  APB.setRandom();
  // make symmetric
  APB = APB.transpose() * APB;
  // make diagonal dominant and A-B positive definite
  Eigen::VectorXd diagonalElements = APB.diagonal();
  diagonalElements = diagonalElements.array().square();
  diagonalElements = diagonalElements.array().sqrt();
  APB *= 1.0e-3;
  APB.diagonal() = diagonalElements;

  // Sigma vector calculators
  auto sigmaVectorCalculator = [&](std::vector<Eigen::MatrixXd>& guessVectors) {
    std::unique_ptr<std::vector<Eigen::MatrixXd>> sigma =
        std::unique_ptr<std::vector<Eigen::MatrixXd>>(new std::vector<Eigen::MatrixXd>(1));
    (*sigma)[0] = APB * guessVectors[0];
    return sigma;
  };

  std::vector<double> frequencies = {1};

  // Build response solver
  auto responseSolver = ResponseSolver(diagonalElements, 1e-6, 100, 100000, frequencies, damping, ppmq, sigmaVectorCalculator);

  // Solutionvectors
  auto itVectors = responseSolver.getEigenvectors();

  // Setup 'exact' lhs-matrix to feed eigen3 with
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(nDim, nDim);
  lhs.block(0 * nDim, 0 * nDim, nDim, nDim) = APB - frequencies[0] * Eigen::MatrixXd::Identity(nDim, nDim);

  // Same for rhs-matrix
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nDim, 3);

  // Determine solution vectors using eigen
  Eigen::MatrixXd solVectors = lhs.colPivHouseholderQr().solve(2 * dipE);

  // Build difference matrix
  Eigen::MatrixXd diffXpY = solVectors - itVectors[0];

  // Compare
  EXPECT_LT(diffXpY.array().cwiseAbs().maxCoeff(), 1e-10);
}

} /* namespace Serenity */
