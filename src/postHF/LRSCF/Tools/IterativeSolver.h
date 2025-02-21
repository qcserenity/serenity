/**
 * @file IterativeSolver.h
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

#ifndef LRSCF_ITERATIVESOLVER
#define LRSCF_ITERATIVESOLVER

/* Include Serenity Internal Headers */
#include "io/FormattedOutput.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace Serenity {

/**
 * @class IterativeSolver IterativeSolver.h
 *
 * @brief Base class for iterative eigenvalue and response solver (not in the current release).
 */
class IterativeSolver {
 public:
  /**
   * @brief Constructor.
   * @param nSets The number of guess vector sets (one and two for Hermitian
   *        and symplectic problems, respectively).
   * @param nDimension Dimension of the response problem.
   * @param nEigen Number of the lowest eigenvalues to be obtained.
   * @param diagonal Orbital-energy difference of the response problem.
   * @param convergenceCriterion If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   * @param maxIterations If not converged after this number of iterations, abort.
   * @param maxSubspaceDimension Will perform a subspace collapse if the number of guess vectors
   *        exceeds this threshold.
   * @param sigmaCalculator A lambda to conveniently form response matrix -- guess vector products.\n
   *        Takes a set of guessvectors as an argument and returns a pointer to the sigmavectors.
   * @param initialGuess The initial guess space might also be passed to the eigenvalue solver.
   */
  IterativeSolver(unsigned nSets, unsigned nDimension, unsigned nEigen, const Eigen::VectorXd& diagonal,
                  double convergenceCriterion, unsigned maxIterations, unsigned maxSubspaceDimension,
                  std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)> sigmaCalculator,
                  std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk,
                  std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess)
    : _nSets(nSets),
      _nDimension(nDimension),
      _nEigen(nEigen),
      _diagonal(diagonal),
      _convergenceCriterion(convergenceCriterion),
      _maxIterations(maxIterations),
      _maxSubspaceDimension(maxSubspaceDimension),
      _sigmaCalculator(sigmaCalculator),
      _writeToDisk(writeToDisk),
      _initialGuess(initialGuess) {
  }

  /**
   * @brief Default destructor.
   */
  virtual ~IterativeSolver() = default;

  /**
   * @brief Triggers the iterative solution procedure.
   */
  void solve() {
    do {
      this->iterate();
      _converged = (_nConverged == _nEigen);
    } while (_converged == false && _nIter < _maxIterations);

    if (_converged) {
      printf("\n  Iterative solver converged in %3i iterations.\n\n", _nIter);
    }
    else if (_maxIterations == 1) {
      printf("\n  Finished FDEc step.\n\n");
    }
    else {
      WarningTracker::printWarning("Warning: Convergence criterion not reached.", true);
    }

    this->postProcessing();
    _nIter = 0;
    _done = true;
  } /* this->solve() */

  /**
   * @brief Returns true if converged.
   */
  bool isConverged() {
    return _converged;
  } /* this->isConverged() */

  /**
   * @brief Returns the current approximate eigenvalues.
   */
  Eigen::VectorXd getEigenvalues() {
    if (!_done) {
      this->solve();
    }
    return _eigenvalues;
  } /* this->getEigenvalues() */

  /**
   * @brief Returns the current approximate eigenvectors.
   */
  std::vector<Eigen::MatrixXd>& getEigenvectors() {
    if (!_done) {
      this->solve();
    }
    return _eigenvectors;
  } /* this->getEigenvectors() */

  /**
   * @brief Performs one iteration.\n
   *        Must be overidden by derived class.
   */
  virtual void iterate() = 0;

  /**
   * @brief Performs operators that are needed after the solution procedure.\n
   *        E.g. obtains X and Y from (X+Y) and (X-Y) and does normalizing.\n
   *        Must be overidden by derived class.
   */
  virtual void postProcessing() = 0;

 protected:
  /**
   * @brief Initializes the eigenvalue solver by printing some info, getting the needed\n
   *        matrices right, calculating the seed and the corresponding sigma vectors.\n
   *        Must be overidden by derived class.
   */
  virtual void initialize() = 0;

  ///@brief Prints some information.
  void printHeader(std::string caption) {
    // Print caption and some general information
    printBigCaption(caption);
    print((std::string) "Number of sets                        : " + _nSets);
    print((std::string) "Dimension of the eigenvalue problem   : " + _nDimension);
    print((std::string) "Number of roots to be determined      : " + _nEigen);
    print((std::string) "Convergence threshold                 : " + _convergenceCriterion);
    print((std::string) "Maximum number of iterations          : " + _maxIterations);
    print((std::string) "Maximum subspace dimension            : " + _maxSubspaceDimension);
  }

  ///@brief Expands the guess space with the current correction vectors.
  void expandSubspace() {
    // Subspace expansion
    unsigned nAppend = 0;
    Eigen::VectorXd norm(_nSets);

    // Orthogonalization against all guess vectors
    for (unsigned iNew = 0; iNew < _correctionVectors[0].cols(); ++iNew) {
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        // oldNorm is the norm of the correction vector before orthogonalization
        double oldNorm = _correctionVectors[iSet].col(iNew).norm();
        for (unsigned j = 0; j < _guessVectors[0].cols(); ++j) {
          double ij = _correctionVectors[iSet].col(iNew).dot(_guessVectors[iSet].col(j));
          double jj = _guessVectors[iSet].col(j).dot(_guessVectors[iSet].col(j));
          _correctionVectors[iSet].col(iNew) -= ij / jj * _guessVectors[iSet].col(j);
        }
        // newNorm is the norm of the part of the correction vector that is orthogonal to all current guess vectors
        double newNorm = _correctionVectors[iSet].col(iNew).norm();
        norm(iSet) = _correctionVectors[iSet].col(iNew).norm() / oldNorm;

        _correctionVectors[iSet].col(iNew) *= (std::max(oldNorm, 1e-7) / newNorm);
      }

      // Determine if to be appended
      if (norm.sum() / std::sqrt(_nSets) > _appendThresh) {
        for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
          _guessVectors[iSet].conservativeResize(_nDimension, _guessVectors[iSet].cols() + 1);
          _guessVectors[iSet].rightCols(1) = _correctionVectors[iSet].col(iNew);
        }
        ++nAppend;
      }
    }

    // Reset correction vectors to only those that got appended
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      _correctionVectors[iSet].conservativeResize(_nDimension, nAppend);
      _correctionVectors[iSet] = _guessVectors[iSet].rightCols(nAppend);
    }

    // Set dimension after expansion
    unsigned newDim = _subDim + nAppend;

    // Perform subspace collapse if needed and calculate new sigma vectors
    if (nAppend > 0 && newDim <= _maxSubspaceDimension) {
      auto newSigma = (*_sigmaCalculator(_correctionVectors));
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        _sigmaVectors[iSet].conservativeResize(_nDimension, newDim);
        _sigmaVectors[iSet].rightCols(nAppend) = newSigma[iSet];
      }
    }
    else {
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        _guessVectors[iSet] = _eigenvectors[iSet];
        _sigmaVectors[iSet] = _sigmaVectors[iSet] * _expansionVectors[iSet];
      }
    }
  } /* this->expandSubspace() */

  ///@brief The number of sets (guess vectors, sigma vectors, eigenvectors, ..).
  unsigned _nSets;

  ///@brief The dimension of the response problem.
  unsigned _nDimension;

  ///@brief The number of roots to be determined.
  unsigned _nEigen;

  ///@brief Orbital-energy difference of the response problem.
  const Eigen::VectorXd& _diagonal;

  /**
   * /@brief If the norm of a residual vector of a root falls beneath this
   *        threshold, the root will be considered converged.
   */
  double _convergenceCriterion;

  ///@brief If not converged after this number of iterations, abort.
  unsigned _maxIterations;

  ///@brief Will perform a subspace collapse if the number of guess vectors exceeds this threshold.
  unsigned _maxSubspaceDimension;

  ///@brief A lambda to conveniently form response matrix -- guess vector products.
  std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)> _sigmaCalculator;

  ///@brief A lambda to store temporary iteration data.
  std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> _writeToDisk;

  ///@brief The initial guess space might also be passed to the eigenvalue solver.
  std::shared_ptr<std::vector<Eigen::MatrixXd>> _initialGuess;

  ///@brief Stores the matrix to be directly diagonalized in each iteration.
  Eigen::MatrixXd _subspaceMatrix;

  ///@brief Stores the subspace eigenvectors from direct diagonalization.
  Eigen::MatrixXd _subspaceEigenvectors;

  ///@brief Stores the subspace eigenvalues from direct diagonalization.
  Eigen::VectorXd _subspaceEigenvalues;

  ///@brief Subspace dimension, equals number of guess/sigma vectors.
  unsigned _subDim;

  ///@brief Subspace metric of guessVectors.
  Eigen::MatrixXd _metric;

  ///@brief Contains _nSets matrices containing the guess vectors.
  std::vector<Eigen::MatrixXd> _guessVectors;

  ///@brief Contains _nSets matrices containing the sigma vectors.
  std::vector<Eigen::MatrixXd> _sigmaVectors;

  ///@brief Contains _nSets matrices containing the current approximate eigenvectors.
  std::vector<Eigen::MatrixXd> _eigenvectors;

  ///@brief Contains _nEigen current approximate eigenvalues.
  Eigen::VectorXd _eigenvalues;

  ///@brief Contains _nSets matrices to store residual vectors for each root.
  std::vector<Eigen::MatrixXd> _residualVectors;

  ///@brief Contains _nSets matrices to store correction vectors obtained from residual vectors.
  std::vector<Eigen::MatrixXd> _correctionVectors;

  ///@brief Contains _nSets matrices to store expansion coefficients obtained from subspace eigenvectors.
  std::vector<Eigen::MatrixXd> _expansionVectors;

  ///@brief The residual norm for each root, which is used to check for convergence.
  Eigen::VectorXd _residualNorms;

  ///@brief Is set to true if the solver is finished, does not necessarily need to be converged.
  bool _done = false;

  ///@brief Is set to true if the solver is converged.
  bool _converged = false;

  /**
   * @brief Correction vectors are orthogonalized against the existing guess space. If their norm after
   *        that is above this threshold, they will be added to the guess space. Is used to prevent linear dependencies.
   */
  double _appendThresh = 1e-3;

  ///@brief Counter for the number of iterations performed.
  unsigned _nIter = 0;

  ///@brief Number of converged roots. If this is equal to _nEigen, the iterative procedure is finished.
  unsigned _nConverged = 0;

  ///@brief Time point at the beginning of an iteration.
  std::chrono::steady_clock::time_point _itStart;

  ///@brief Time point at the end of an iteration.
  std::chrono::steady_clock::time_point _itEnd;

}; /* class IterativeSolver */
} /* namespace Serenity */

#endif /* LRSCF_ITERATIVESOLVER */
