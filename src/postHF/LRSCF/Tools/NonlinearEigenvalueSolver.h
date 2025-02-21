/**
 * @file NonlinearEigenvalueSolver.h
 *
 * @date March 18, 2020
 * @author Niklas Niemeyer
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

#ifndef NONLINEAREIGENVALUESOLVER
#define NONLINEAREIGENVALUESOLVER

/* Include Serenity Internal Headers */
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <chrono>
#include <memory>

namespace Serenity {
/**
 * @class NonlinearEigenvalueSolver
 * @brief A class for the iterative solution of the partitioned response eigenvalue problem.\n
 *
 * This class includes solution techniques for Hermitian and non-Hermitian eigenvalue\n
 * problems, as required by ADC(2), CIS(D_inf) and CC2.
 */
class NonlinearEigenvalueSolver {
 public:
  /**
   * @brief Constructor.
   * @param nEigen Number of eigenvalues to be determined.
   * @param diagonal Orbital-energy differences (zeroth-order Jacobian).
   * @param conv Convergence threshold for residuals.
   * @param preopt Preoptimization convergence threshold for residuals.
   * @param diis DIIS to be used after preoptimization or not.
   * @param diisStore Maximum number of DIIS vectors to be stored.
   * @param maxCycles Maximum number of cycles to be performed in each procedure.
   * @param rootFollowing Tells the eigenvalue solver to follow the input or not.
   * @param method Excited-state wavefunction method (CC2, CIS(Dinf), CIS(D) or ADC(2)).
   * @param eigenvectors Eigenvectors to be optimized (passed by reference).
   * @param eigenvalues Eigenvalues to be optimized (passed by reference).
   * @param sigmaCalculator Lambda function to perform matrix--vector multiplication.
   * @param writeToDisk Lambda function to store unconverged solutions on disk.
   */
  NonlinearEigenvalueSolver(
      unsigned nEigen, const Eigen::VectorXd& diagonal, double conv, double preopt, bool diis, unsigned diisStore,
      unsigned maxCycles, bool rootFollowing, Options::LR_METHOD method, Eigen::MatrixXd& eigenvectors,
      Eigen::VectorXd& eigenvalues,
      std::function<std::unique_ptr<Eigen::MatrixXd>(Eigen::MatrixXd& guessVectors, Eigen::VectorXd& guessValues)> sigmaCalculator,
      std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk = [](std::vector<Eigen::MatrixXd>&,
                                                                                            Eigen::VectorXd&) {});

  /**
   * @brief Default destructor
   */
  virtual ~NonlinearEigenvalueSolver() = default;

  /**
   * @brief Triggers the macro-solution procedure of the non-linear eigenvalue solver.
   */
  void solve();

 private:
  /**
   * @brief Performs one iteration of the preoptimization.
   */
  void iteratePreopt();

  /**
   * @brief Performs one iteration of the final DIIS cycles.
   */
  void iterateDIIS();

  /**
   * @brief Performs the subspace expansion.
   */
  void expandSubspace();

  /**
   * @brief Really only prints the currently stored solution.
   */
  void postProcessing();

  ///@brief The number of roots to be determined.
  unsigned _nEigen;

  ///@brief Orbital-energy difference of the response problem.
  const Eigen::VectorXd& _diagonal;

  /**
   * @brief If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   */
  double _conv;

  ///@brief If DIIS used, preopt convergence threshold.
  double _preopt;

  //@brief Determines if a DIIS is to be used after the preoptimization procedure.
  bool _diis;

  ///@brief Restricts the DIIS vectors to be stored.
  unsigned _diisStore;

  ///@brief If not converged after this number of iterations, abort.
  unsigned _maxCycles;

  ///@brief Tells the eigenvalue solver to follow in the initial eigenvectors or not.
  bool _rootFollowing;

  ///@brief The underlying excited-state method.
  Options::LR_METHOD _method;

  ///@brief Contains _nSets matrices containing the current approximate eigenvectors.
  Eigen::MatrixXd& _eigenvectors;

  ///@brief Contains _nEigen current approximate eigenvalues.
  Eigen::VectorXd& _eigenvalues;

  ///@brief A lambda to conveniently form response matrix -- guess vector products.
  std::function<std::unique_ptr<Eigen::MatrixXd>(Eigen::MatrixXd&, Eigen::VectorXd&)> _sigmaCalculator;

  ///@brief A lambda to store temporary iteration data.
  std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> _writeToDisk;

  ///@brief The dimension of the response problem.
  unsigned _nDimension;

  ///@brief Will perform a subspace collapse if the number of guess vectors exceeds this threshold.
  unsigned _maxSubspaceDimension;

  ///@brief Is set to true if the solver is finished, does not necessarily need to be converged.
  bool _done;

  ///@brief Convergence threshold for the current micro iteration.
  double _microConv;

  ///@brief Is set to true if the solver is converged.
  bool _converged;

  ///@brief Bool indicating whether a new conv treshold needs to be applied.
  bool _macroConv;

  ///@brief Counter for the number of iterations performed.
  unsigned _nIter;

  ///@brief Number of converged roots. If this is equal to _nEigen, the iterative procedure is finished.
  unsigned _nConverged;

  ///@brief The current best eigenvalues
  Eigen::VectorXd _newEigenvalues;

  ///@brief Used to store the eigenvalues for which sigmavectors are to be calculated.
  Eigen::VectorXd _sigmaEigenvalues;

  ///@brief Contains _nSets matrices containing the guess vectors.
  Eigen::MatrixXd _guessVectors;

  ///@brief Contains the initial guess vectors.
  Eigen::MatrixXd _initialGuess;

  ///@brief Contains _nSets matrices containing the sigma vectors.
  Eigen::MatrixXd _sigmaVectors;

  ///@brief Contains _nSets matrices to store residual vectors for each root.
  Eigen::MatrixXd _residualVectors;

  ///@brief Contains _nSets matrices to store correction vectors obtained from residual vectors.
  Eigen::MatrixXd _correctionVectors;

  ///@brief Contains _nSets matrices to store expansion coefficients obtained from subspace eigenvectors.
  Eigen::MatrixXd _expansionVectors;

  ///@brief The residual norm for each root, which is used to check for convergence.
  Eigen::VectorXd _residualNorms;

  ///@brief Storage of amplitudes for DIIS for each root.
  std::vector<Eigen::MatrixXd> _amplitudeStorage;

  ///@brief Storage of residuals for DIIS for each root.
  std::vector<Eigen::MatrixXd> _residualStorage;

  ///@brief Time point at the beginning of an iteration.
  std::chrono::steady_clock::time_point _itStart;

  ///@brief Time point at the end of an iteration.
  std::chrono::steady_clock::time_point _itEnd;

}; /* class NonlinearEigenvalueSolver */
} /* namespace Serenity */

#endif /* NONLINEAREIGENVALUESOLVER */