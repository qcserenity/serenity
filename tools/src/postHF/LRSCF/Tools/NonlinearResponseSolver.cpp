/**
 * @file NonlinearResponseSolver.cpp
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

/* Include Class Header*/
#include "postHF/LRSCF/Tools/NonlinearResponseSolver.h"
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "parameters/Constants.h"
#include "settings/LRSCFOptions.h"

namespace Serenity {

NonlinearResponseSolver::NonlinearResponseSolver(const Eigen::VectorXd& diagonal, double convergenceCriterion,
                                                 unsigned maxIterations, unsigned maxSubspaceDimension,
                                                 std::vector<double> frequencies, std::vector<Eigen::MatrixXd> ppmq,
                                                 NonlinearSigmaCalculator sigmaCalculator)
  : _diagonal(diagonal),
    _conv(convergenceCriterion),
    _maxCycles(maxIterations),
    _nDimension(diagonal.size()),
    _maxSubspaceDimension(maxSubspaceDimension),
    _frequencies(frequencies),
    _nFreqs(_frequencies.size()),
    _ppmq(ppmq),
    _sigmaCalculator(sigmaCalculator) {
  this->initialize();
}

void NonlinearResponseSolver::initialize() {
  printBigCaption("Nonlinear Response Solver");

  _nEigen = 0;
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    _nEigen += _ppmq[iFreq].cols();
  }

  printf("    Number of roots to be determined  : %-5i\n", _nEigen);
  printf("    Maximum number of iterations      : %-5i\n", _maxCycles);
  printf("    Maximum subspace dimension        : %-5i\n", _maxSubspaceDimension);
  printf("    Size of single-excitation space   : %-5i\n", _nDimension);
  printf("    Convergence threshold             : %-5.1e\n", _conv);

  // Set dimensions
  _guessVectors.resize(_nFreqs);
  _sigmaVectors.resize(_nFreqs);
  _eigenvectors.resize(_nFreqs);
  _residualVectors.resize(_nFreqs);
  _correctionVectors.resize(_nFreqs);
  _residualNorms.resize(_nFreqs);
  _freqConverged.resize(_nFreqs);

  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    _guessVectors[iFreq] = Eigen::MatrixXd::Zero(_nDimension, _ppmq[iFreq].cols());
    _eigenvectors[iFreq] = Eigen::MatrixXd::Zero(_nDimension, _ppmq[iFreq].cols());
    _residualVectors[iFreq] = Eigen::MatrixXd::Zero(_nDimension, _ppmq[iFreq].cols());
    _correctionVectors[iFreq] = Eigen::MatrixXd::Zero(_nDimension, _ppmq[iFreq].cols());
    _freqConverged[iFreq] = false;
  }

  this->seed();

  Orthogonalization::modifiedGramSchmidtLinDep(_guessVectors, 1e-4);

  _itStart = std::chrono::steady_clock::now();
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    _sigmaVectors[iFreq] =
        (*_sigmaCalculator(_guessVectors[iFreq], Eigen::VectorXd::Constant(_ppmq[iFreq].cols(), _frequencies[iFreq])));
  }
}

void NonlinearResponseSolver::solve() {
  _nIter = 1;
  do {
    this->iterate();
    ++_nIter;
  } while (_converged == false && _nIter <= _maxCycles);

  if (_converged) {
    _done = true;
    printf("\n  Nonlinear response converged in %3i iterations.\n\n", _nIter - 1);
  }
  else {
    WarningTracker::printWarning("Warning: Convergence criterion not reached.", true);
  }
}

void NonlinearResponseSolver::iterate() {
  _nConverged = 0;
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    // Check if this frequency is already finished.
    if (_freqConverged[iFreq]) {
      _nConverged += _ppmq[iFreq].cols();
      continue;
    }

    // Construct subspace problem.
    Eigen::MatrixXd lhs =
        _guessVectors[iFreq].transpose() * (_sigmaVectors[iFreq] - _frequencies[iFreq] * _guessVectors[iFreq]);
    Eigen::MatrixXd rhs = _guessVectors[iFreq].transpose() * _ppmq[iFreq];

    // Ritz vectors.
    Eigen::MatrixXd expansionVectors = lhs.fullPivLu().solve(rhs);
    _eigenvectors[iFreq] = _guessVectors[iFreq] * expansionVectors;

    // Residual vectors.
    _residualVectors[iFreq] =
        (_sigmaVectors[iFreq] * expansionVectors) - (_frequencies[iFreq] * _eigenvectors[iFreq]) - _ppmq[iFreq];
    _residualNorms[iFreq] = _residualVectors[iFreq].colwise().norm();

    // Convergence check.
    unsigned counter = 0;
    for (unsigned i = 0; i < _ppmq[iFreq].cols(); ++i) {
      if (_residualNorms[iFreq](i) < _conv) {
        ++counter;
        ++_nConverged;
      }
    }
    if (counter == _ppmq[iFreq].cols()) {
      _freqConverged[iFreq] = true;
    }
  }

  _itEnd = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();

  Eigen::VectorXd maxNorms(_nFreqs);
  int subDim = 0;
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    maxNorms(iFreq) = _residualNorms[iFreq].maxCoeff();
    subDim += _guessVectors[iFreq].cols();
  }
  unsigned iMax;
  maxNorms.maxCoeff(&iMax);
  // Print Iteration Data
  if (_nIter == 1) {
    printTableHead("\n   it.    dimension    time (min)    converged     max. norm");
  }
  printf("%5i %10i %14.3f %11i %13.2e  (%2i)\n", _nIter, subDim, duration / 60.0, _nConverged, maxNorms.maxCoeff(), iMax + 1);

  _itStart = std::chrono::steady_clock::now();

  // Skip the rest of this iteration if all solution vectors converged
  _converged = (_nConverged == _nEigen);
  if (_converged || _nIter > _maxCycles) {
    return;
  }

  this->precondition();
  this->expandSubspace();
}

void NonlinearResponseSolver::seed() {
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    for (unsigned i = 0; i < _ppmq[iFreq].cols(); ++i) {
      for (unsigned ia = 0; ia < _nDimension; ++ia) {
        _guessVectors[iFreq](ia, i) = 1 / (_frequencies[iFreq] - _diagonal(ia)) * _ppmq[iFreq](ia, i);
      }
    }
  }
}

void NonlinearResponseSolver::precondition() {
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    _correctionVectors[iFreq].resize(_nDimension, 0);
    unsigned index = 0;
    for (unsigned i = 0; i < _ppmq[iFreq].cols(); ++i) {
      if (_residualNorms[iFreq](i) < _conv) {
        continue;
      }
      _correctionVectors[iFreq].conservativeResize(Eigen::NoChange, _correctionVectors[iFreq].cols() + 1);
      for (unsigned ia = 0; ia < _nDimension; ++ia) {
        _correctionVectors[iFreq](ia, index) = 1 / (_frequencies[iFreq] - _diagonal(ia)) * _residualVectors[iFreq](ia, i);
      }
      ++index;
    }
  }
}

void NonlinearResponseSolver::expandSubspace() {
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    unsigned nAppend = 0;
    for (unsigned iNew = 0; iNew < _correctionVectors[iFreq].cols(); ++iNew) {
      _correctionVectors[iFreq].col(iNew).normalize();
      for (unsigned iGuess = 0; iGuess < _guessVectors[iFreq].cols(); ++iGuess) {
        double ij = _correctionVectors[iFreq].col(iNew).dot(_guessVectors[iFreq].col(iGuess));
        double jj = _guessVectors[iFreq].col(iGuess).dot(_guessVectors[iFreq].col(iGuess));
        _correctionVectors[iFreq].col(iNew) -= ij / jj * _guessVectors[iFreq].col(iGuess);
      }
      _correctionVectors[iFreq].col(iNew).normalize();
      _guessVectors[iFreq].conservativeResize(Eigen::NoChange, _guessVectors[iFreq].cols() + 1);
      _guessVectors[iFreq].rightCols(1) = _correctionVectors[iFreq].col(iNew);
      ++nAppend;
    }
    _correctionVectors[iFreq] = _guessVectors[iFreq].rightCols(nAppend);

    // New sigma vectors.
    if (nAppend > 0) {
      Eigen::MatrixXd newSigma =
          (*_sigmaCalculator(_correctionVectors[iFreq], Eigen::VectorXd::Constant(nAppend, _frequencies[iFreq])));
      _sigmaVectors[iFreq].conservativeResize(Eigen::NoChange, _sigmaVectors[iFreq].cols() + nAppend);
      _sigmaVectors[iFreq].rightCols(nAppend) = newSigma;
    }
  }
}

std::vector<Eigen::MatrixXd> NonlinearResponseSolver::getEigenvectors() {
  if (_done == false) {
    this->solve();
  }
  return _eigenvectors;
}

} /* namespace Serenity */
