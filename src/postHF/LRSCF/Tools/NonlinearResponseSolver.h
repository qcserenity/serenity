/**
 * @file NonlinearResponseSolver.h
 *
 * @date May 30, 2019
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

#ifndef LRSCF_NONLINEARRESPONSESOLVER
#define LRSCF_NONLINEARRESPONSESOLVER

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/IterativeSolver.h"
#include "postHF/LRSCF/Tools/SigmaCalculator.h"

namespace Serenity {

/**
 * @class NonlinearResponseSolver
 * @brief This class is the analog of the NonlinearEigenvalueSolver, which diagonalizes a non-linear eigenvalue
 * problem (the matrix is a function of its eigenvalue). This class, while solving a linear system, is non-linear
 * in the sense that the matrix to form-matrix-vector products with depends on the frequency employed (as required
 * in the response context).
 */
class NonlinearResponseSolver {
 public:
  /**
   * @brief Constructor.
   * @param diagonal Orbital-energy difference of the CC2 problem.
   * @param convergenceCriterion If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   * @param maxIterations If not converged after this number of iterations, abort.
   * @param maxSubspaceDimension Will perform a subspace collapse if the number of guess vectors
   *        exceeds this threshold.
   * @param frequencies The frequencies to be solved for.
   * @param ppmq The right-hand side of the linear system.
   * @param sigmaCalculator A lambda to conveniently form response matrix -- guess vector products.\n
   *        Takes a set of guessvectors as an argument and returns a pointer to the sigmavectors.
   */
  NonlinearResponseSolver(const Eigen::VectorXd& diagonal, double convergenceCriterion, unsigned maxIterations,
                          unsigned maxSubspaceDimension, std::vector<double> frequencies,
                          std::vector<Eigen::MatrixXd> ppmq, NonlinearSigmaCalculator sigmaCalculator);

  /**
   * @brief Default destructor
   */
  virtual ~NonlinearResponseSolver() = default;

  /**
   * @brief Triggers the solution procedure of the non-linear response solver.
   */
  void solve();

  /**
   * @brief Returns the solution vectors.
   */
  std::vector<Eigen::MatrixXd> getEigenvectors();

 private:
  /**
   * @brief Initialized the solver.
   */
  void initialize();

  /**
   * @brief Creates initial guess vectors.
   */
  void seed();

  /**
   * @brief Performs one iteration.
   */
  void iterate();

  /**
   * @brief Preconditions correction vectors with diagonal approximation.
   */
  void precondition();

  /**
   * @brief Performs the subspace expansion.
   */
  void expandSubspace();

  ///@brief Orbital-energy difference of the response problem.
  const Eigen::VectorXd& _diagonal;

  /**
   * @brief If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   */
  double _conv;

  ///@brief If not converged after this number of iterations, abort.
  unsigned _maxCycles;

  ///@brief Contains approximate solution vectors.
  std::vector<Eigen::MatrixXd> _eigenvectors;

  ///@brief The dimension of the response problem.
  unsigned _nDimension;

  ///@brief Will perform a subspace collapse if the number of guess vectors exceeds this threshold.
  unsigned _maxSubspaceDimension;

  ///@brief The frequencies to be solved for.
  std::vector<double> _frequencies;

  ///@brief The number of frequencies considered.
  unsigned _nFreqs;

  ///@brief The number of roots to be determined.
  unsigned _nEigen;

  ///@brief The right-hand sides to be solved for.
  std::vector<Eigen::MatrixXd> _ppmq;

  ///@brief A lambda to conveniently form response matrix -- guess vector products.
  NonlinearSigmaCalculator _sigmaCalculator;

  ///@brief Is set to true if the solver is finished, does not necessarily need to be converged.
  bool _done = false;

  ///@brief Is set to true if the solver is converged.
  bool _converged;

  ///@brief Counter for the number of iterations performed.
  unsigned _nIter;

  ///@brief Number of converged roots. If this is equal to _nEigen, the iterative procedure is finished.
  unsigned _nConverged;

  ///@brief Indicates whether this frequency is converged.
  std::vector<bool> _freqConverged;

  ///@brief Contains _nSets matrices containing the guess vectors.
  std::vector<Eigen::MatrixXd> _guessVectors;

  ///@brief Contains _nSets matrices containing the sigma vectors.
  std::vector<Eigen::MatrixXd> _sigmaVectors;

  ///@brief Contains _nSets matrices to store residual vectors for each root.
  std::vector<Eigen::MatrixXd> _residualVectors;

  ///@brief Contains _nFreqs matrices to store correction vectors obtained from residual vectors.
  std::vector<Eigen::MatrixXd> _correctionVectors;

  ///@brief The residual norm for each root, which is used to check for convergence.
  std::vector<Eigen::VectorXd> _residualNorms;

  ///@brief Time point at the beginning of an iteration.
  std::chrono::steady_clock::time_point _itStart;

  ///@brief Time point at the end of an iteration.
  std::chrono::steady_clock::time_point _itEnd;

}; /* class NonlinearResponseSolver */
} /* namespace Serenity */

#endif /* LRSCF_NONLINEARRESPONSESOLVER */
