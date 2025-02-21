/**
 * @file NonlinearEigenvalueSolver.cpp
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

/* Include Class Header*/
#include "postHF/LRSCF/Tools/NonlinearEigenvalueSolver.h"
/* Include Serenity Internal Headers */
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"

namespace Serenity {

NonlinearEigenvalueSolver::NonlinearEigenvalueSolver(
    unsigned nEigen, const Eigen::VectorXd& diagonal, double conv, double preopt, bool diis, unsigned diisStore,
    unsigned maxCycles, bool rootFollowing, Options::LR_METHOD method, Eigen::MatrixXd& eigenvectors,
    Eigen::VectorXd& eigenvalues,
    std::function<std::unique_ptr<Eigen::MatrixXd>(Eigen::MatrixXd& guessVectors, Eigen::VectorXd& guessValues)> sigmaCalculator,
    std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk)
  : _nEigen(nEigen),
    _diagonal(diagonal),
    _conv(conv),
    _preopt(preopt),
    _diis(diis),
    _diisStore(diisStore),
    _maxCycles(maxCycles),
    _rootFollowing(rootFollowing),
    _method(method),
    _eigenvectors(eigenvectors),
    _eigenvalues(eigenvalues),
    _sigmaCalculator(sigmaCalculator),
    _writeToDisk(writeToDisk),
    _nDimension(diagonal.size()),
    _maxSubspaceDimension(diisStore * _nEigen),
    _done(false),
    _converged(false),
    _macroConv(true),
    _nIter(0),
    _nConverged(0) {
  printBigCaption("Nonlinear Eigenvalue Solver");

  // If diis to be used, initialize the diis matrices for each root.
  if (_diis) {
    _amplitudeStorage = std::vector<Eigen::MatrixXd>(_nEigen, Eigen::MatrixXd::Zero(0, 0));
    _residualStorage = std::vector<Eigen::MatrixXd>(_nEigen, Eigen::MatrixXd::Zero(0, 0));
  }
  // Otherwise, make the quasi-linear davidson converge up to conv.
  else {
    _preopt = _conv;
  }

  // print information.
  printf("    Number of roots to be determined  : %-5i\n", _nEigen);
  printf("    Maximum number of iterations      : %-5i\n", _maxCycles);
  printf("    Maximum subspace dimension        : %-5i\n", _maxSubspaceDimension);
  printf("    Size of single-excitation space   : %-5i\n", _nDimension);
  if (_diis) {
    printf("    Preopt Davidson threshold         : %-5.1e\n", _preopt);
    printf("    DIIS convergence threshold        : %-5.1e\n\n", _conv);
  }
  else {
    printf("    Convergence threshold             : %-5.1e\n\n", _conv);
  }

  _initialGuess = _eigenvectors;
  _guessVectors = _eigenvectors;

  _itStart = std::chrono::steady_clock::now();
  _sigmaVectors = *_sigmaCalculator(_guessVectors, _eigenvalues);
}

void NonlinearEigenvalueSolver::solve() {
  // CIS(D) is a one-shot procedure, so sigmavectors only have to be calculated once.
  if (_method == Options::LR_METHOD::CISD) {
    _eigenvalues = (_eigenvectors.transpose() * _sigmaVectors).diagonal();
    _residualNorms = (_sigmaVectors - _eigenvectors * _eigenvalues.asDiagonal()).colwise().norm();

    _itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();
    unsigned iMax;
    double maxResNorm = _residualNorms.maxCoeff(&iMax);
    printf("   it.   dim.   time (min)    conv.    conv thr        max norm     \n");
    printf("  ------------------------------------------------------------------\n");
    printf("%5i %6i %11.3f %8i %13.2e %12.2e   (%2i)\n", ++_nIter, (int)_guessVectors.cols(), duration / 60.0, _nEigen,
           0.0, maxResNorm, iMax + 1);
  }
  /**
   * Otherwise: first preoptimize roots with a quasi-linear davidson up to a convergence
   * threshold of preopt, and then use a diis extrapolation for full convergence to conv.
   *
   * If diis is turned off, simply converge to conv using the Davidson algorithm.
   */
  else {
    _nIter = 0;
    do {
      this->iteratePreopt();
      ++_nIter;
      if (_macroConv) {
        auto disk = std::vector<Eigen::MatrixXd>(1, _eigenvectors);
        _writeToDisk(disk, _eigenvalues);
      }
    } while (_converged == false && _nIter < _maxCycles);

    if (_maxCycles == 1) {
      printf("\n  Finished FDEc step.\n\n");
    }
    else if (_converged) {
      printf("\n  Nonlinear eigensolver converged in %3i iterations.\n\n", _nIter);
    }
    else {
      WarningTracker::printWarning("Warning: Convergence criterion not reached.", true);
    }
    if (_diis) {
      _converged = false;
      printf("\n  Consider these results preoptimized!\n");
      this->postProcessing();

      printBigCaption("DIIS Eigenvalue Solver");

      _nIter = 0;
      do {
        this->iterateDIIS();
        ++_nIter;
      } while (_converged == false && _nIter < _maxCycles);
      if (_converged) {
        printf("\n  DIIS Eigensolver converged in %3i iterations.\n\n", _nIter);
      }
      else {
        WarningTracker::printWarning("Warning: Convergence criterion not reached.", true);
      }
    }
  }

  this->postProcessing();
}

void NonlinearEigenvalueSolver::iteratePreopt() {
  // Solve subspace problem.
  Eigen::MatrixXd subspaceMatrix = _guessVectors.transpose() * _sigmaVectors;
  Eigen::MatrixXd metricMatrix = _guessVectors.transpose() * _guessVectors;

  Eigen::EigenSolver<Eigen::MatrixXd> subspaceSolver(metricMatrix.inverse() * subspaceMatrix);
  Eigen::MatrixXd subspaceEigenvectors = subspaceSolver.eigenvectors().real();
  Eigen::VectorXd subspaceEigenvalues = subspaceSolver.eigenvalues().real();

  // Sort eigenvalues in ascending order.
  unsigned iMin;
  for (unsigned i = 0; i < subspaceEigenvalues.size(); ++i) {
    subspaceEigenvalues.tail(subspaceEigenvalues.size() - i).minCoeff(&iMin);
    subspaceEigenvalues.row(i).swap(subspaceEigenvalues.row(iMin + i));
    subspaceEigenvectors.col(i).swap(subspaceEigenvectors.col(iMin + i));
  }

  _newEigenvalues = subspaceEigenvalues;
  _expansionVectors = subspaceEigenvectors;

  unsigned iMax;
  // If only one iteration is requested, follow the input roots to preserve the ordering.
  // After that, it's good to just have the eigenpairs sorted in ascending order.
  if (_rootFollowing) {
    Eigen::MatrixXd eigenvectors = _guessVectors * _expansionVectors;

    // Root following: follow those roots resembling the initial guess best.
    Eigen::MatrixXd overlap = (eigenvectors.transpose() * _initialGuess).cwiseAbs();

    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      overlap.col(iEigen).maxCoeff(&iMax);
      _newEigenvalues.row(iEigen).swap(_newEigenvalues.row(iMax));
      _expansionVectors.col(iEigen).swap(_expansionVectors.col(iMax));
    }
  }

  // Only take the relevant roots.
  _newEigenvalues = _newEigenvalues.head(_nEigen).eval();
  _expansionVectors = _expansionVectors.leftCols(_nEigen).eval();

  double maxEnergyDiff = (_eigenvalues - _newEigenvalues).cwiseAbs().maxCoeff();

  if (_macroConv) {
    _microConv = std::max(maxEnergyDiff, _preopt);
    _macroConv = false;
  }
  else if (maxEnergyDiff > _microConv) {
    _microConv = std::max(maxEnergyDiff, _preopt);
  }

  // Ritz vectors.
  _eigenvectors = _guessVectors * _expansionVectors;

  // Residual vectors and norms.
  _residualVectors = _sigmaVectors * _expansionVectors - _eigenvectors * _newEigenvalues.asDiagonal();
  _residualNorms = _residualVectors.colwise().norm();

  // Check preopt convergence.
  _nConverged = 0;
  for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
    if (_residualNorms(iEigen) < _microConv) {
      ++_nConverged;
    }
  }

  _itEnd = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();

  double maxResNorm = _residualNorms.maxCoeff(&iMax);
  if (_nIter == 0) {
    printf("   it.   dim.   time (min)    conv.    conv thr        max norm     \n");
    printf("  ------------------------------------------------------------------\n");
  }
  printf("%5i %6i %11.3f %8i %13.2e %12.2e   (%2i)\n", _nIter + 1, (int)_guessVectors.cols(), duration / 60.0,
         _nConverged, _microConv, maxResNorm, iMax + 1);

  // This is for the FDEc case.
  if (_maxCycles == 1) {
    printf("\n");
    printSmallCaption("Subspace Matrix / eV");
    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      printf(" ");
      for (unsigned jEigen = 0; jEigen < _nEigen; ++jEigen) {
        printf("%12.3e", subspaceMatrix(iEigen, jEigen) * HARTREE_TO_EV);
      }
      printf("\n");
    }
    printf("\n");
    printSmallCaption("Metric Matrix / eV");
    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      printf(" ");
      for (unsigned jEigen = 0; jEigen < _nEigen; ++jEigen) {
        printf("%12.3e", metricMatrix(iEigen, jEigen));
      }
      printf("\n");
    }
    printf("\n");

    _eigenvectors.colwise().normalize();
    _eigenvalues = _newEigenvalues;
    return;
  }

  _itStart = std::chrono::steady_clock::now();

  // Evaluate convergence.
  _macroConv = (_nConverged == _nEigen);
  if (_macroConv) {
    if (_nIter == 0) {
      _eigenvectors = _guessVectors;
    }
    _eigenvectors.colwise().normalize();
    _eigenvalues = _newEigenvalues;

    if (maxResNorm < _preopt && std::abs(_microConv - _preopt) < 1e-12) {
      _converged = true;
    }
    else {
      _guessVectors = _eigenvectors;
      _sigmaVectors = *_sigmaCalculator(_guessVectors, _eigenvalues);
      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    }
  }
  else {
    // Reset correction vectors.
    _correctionVectors = Eigen::MatrixXd::Zero(_nDimension, _nEigen - _nConverged);
    _sigmaEigenvalues = Eigen::VectorXd::Zero(_nEigen - _nConverged);

    // Use diagonal approximation to precondition residual vectors.
    unsigned iCount = 0;
    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      if (_residualNorms(iEigen) < _microConv) {
        continue;
      }
      _sigmaEigenvalues(iCount) = _eigenvalues(iEigen);
      for (unsigned ia = 0; ia < _nDimension; ++ia) {
        _correctionVectors(ia, iCount) -= _residualVectors(ia, iEigen) / (_diagonal(ia) - _eigenvalues(iEigen));
      }
      ++iCount;
    }

    this->expandSubspace();
  }
} /* this->iteratePreopt() */

void NonlinearEigenvalueSolver::iterateDIIS() {
  // Get sigmavectors for all roots.
  if (_nIter == 0) {
    _eigenvectors.colwise().normalize();
    _sigmaVectors.resize(_nDimension, _nEigen);
    _correctionVectors = _eigenvectors;
    _sigmaEigenvalues = _eigenvalues;

    _itStart = std::chrono::steady_clock::now();
  }

  // Only calculate sigmavectors for unconverged roots.
  Eigen::MatrixXd newSigma = *_sigmaCalculator(_correctionVectors, _sigmaEigenvalues);

  // Distribute newly calculated sigma vectors.
  unsigned iCount = 0;
  for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
    if (_residualNorms(iEigen) < _conv) {
      continue;
    }
    _sigmaVectors.col(iEigen) = newSigma.col(iCount);
    ++iCount;
  }

  // Just store difference in new.
  _newEigenvalues = _eigenvalues - (_eigenvectors.transpose() * _sigmaVectors).diagonal();

  // Get new eigenvalues.
  _eigenvalues = (_eigenvectors.transpose() * _sigmaVectors).diagonal();
  _residualVectors = (_sigmaVectors - _eigenvectors * _eigenvalues.asDiagonal());
  _residualNorms = _residualVectors.colwise().norm();

  // Evaluate convergence.
  _nConverged = 0;
  for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
    if (_residualNorms(iEigen) < _conv) {
      ++_nConverged;
    }
  }

  // Return if converged.
  if (_nConverged == _nEigen) {
    _converged = true;
  }
  // Repeat diis for unconverged roots.
  else {
    _correctionVectors = Eigen::MatrixXd::Zero(_nDimension, _nEigen - _nConverged);
    _sigmaEigenvalues = Eigen::VectorXd::Zero(_nEigen - _nConverged);

    iCount = 0;
    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      if (_residualNorms(iEigen) < _conv) {
        continue;
      }
      for (unsigned ia = 0; ia < _nDimension; ++ia) {
        _eigenvectors(ia, iEigen) -= _residualVectors(ia, iEigen) / (_diagonal(ia) - _eigenvalues(iEigen));
      }
      // DIIS.
      if (_nIter < _diisStore) {
        _residualStorage[iEigen].conservativeResize(_nDimension, _residualStorage[iEigen].cols() + 1);
        _amplitudeStorage[iEigen].conservativeResize(_nDimension, _amplitudeStorage[iEigen].cols() + 1);
      }
      _residualStorage[iEigen].col(_nIter % _diisStore) = _residualVectors.col(iEigen);
      _amplitudeStorage[iEigen].col(_nIter % _diisStore) = _eigenvectors.col(iEigen);

      Eigen::VectorXd rhs(Eigen::VectorXd::Zero(_nIter + 2));
      Eigen::MatrixXd lhs(Eigen::MatrixXd::Zero(_nIter + 2, _nIter + 2));

      for (unsigned iStored = 0; iStored <= _nIter; ++iStored) {
        lhs(iStored, _nIter + 1) = -1.0;
        lhs(_nIter + 1, iStored) = -1.0;
      }
      lhs.topLeftCorner(_nIter + 1, _nIter + 1) = _residualStorage[iEigen].transpose() * _residualStorage[iEigen];
      rhs(_nIter + 1) = -1.0;
      _eigenvectors.col(iEigen) = _amplitudeStorage[iEigen] * lhs.householderQr().solve(rhs).head(_nIter + 1);
      _eigenvectors.col(iEigen).normalize();

      // Build vector and energy for sigmavector for this unconverged root.
      _sigmaEigenvalues(iCount) = _eigenvalues(iEigen);
      _correctionVectors.col(iCount) = _eigenvectors.col(iEigen);
      ++iCount;
    } /* iEigen */
  }

  // Print iteration data.
  _itEnd = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();

  unsigned iMax;
  _residualNorms.maxCoeff(&iMax);
  if (_nIter == 0) {
    printf("   it.   time (min)    conv.    max energy diff      max norm    \n");
    printf("  ---------------------------------------------------------------\n");
  }
  printf("%5i %11.3f %8i %17.2e %13.2e   (%2i)\n", _nIter + 1, duration / 60.0, _nConverged,
         _newEigenvalues.cwiseAbs().maxCoeff(), _residualNorms.maxCoeff(), iMax + 1);

  _itStart = std::chrono::steady_clock::now();
} /* this->iterateDIIS() */

void NonlinearEigenvalueSolver::expandSubspace() {
  unsigned nAppend = 0;

  // Normalize.
  for (unsigned iNew = 0; iNew < _correctionVectors.cols(); ++iNew) {
    _guessVectors.conservativeResize(_nDimension, _guessVectors.cols() + 1);
    _guessVectors.rightCols(1) = _correctionVectors.col(iNew) / _correctionVectors.col(iNew).norm();
    ++nAppend;
  }

  // Reset correction vectors to only those that got appended.
  _correctionVectors = _guessVectors.rightCols(nAppend);

  // Set dimension after expansion.
  unsigned newDim = _guessVectors.cols();

  // Calculate new sigma vectors.
  Eigen::MatrixXd newSigma = (*_sigmaCalculator(_correctionVectors, _sigmaEigenvalues));
  // If there is more space, then just append new sigma vectors.
  if (newDim <= _maxSubspaceDimension) {
    _sigmaVectors.conservativeResize(_nDimension, newDim);
    _sigmaVectors.rightCols(nAppend) = newSigma;
  }
  // Otherwise, collapse subspace and then append new sigma vectors.
  else {
    // Set guessVectors to current Ritz vectors.
    _guessVectors = _eigenvectors;
    _sigmaVectors = _sigmaVectors * _expansionVectors;

    // Don't neglect current correction vectors.
    _guessVectors.conservativeResize(_nDimension, _nEigen + nAppend);
    _guessVectors.rightCols(nAppend) = _correctionVectors;

    // Append corresponding sigma vectors.
    _sigmaVectors.conservativeResize(_nDimension, _nEigen + nAppend);
    _sigmaVectors.rightCols(nAppend) = newSigma;
  }
}

void NonlinearEigenvalueSolver::postProcessing() {
  if (_method != Options::LR_METHOD::CISD) {
    // Sort eigenpairs in ascending order one last time.
    unsigned iMin;
    for (unsigned i = 0; i < _eigenvalues.size(); ++i) {
      _eigenvalues.tail(_eigenvalues.size() - i).minCoeff(&iMin);
      _eigenvalues.row(i).swap(_eigenvalues.row(iMin + i));
      _eigenvectors.col(i).swap(_eigenvectors.col(iMin + i));
      _residualNorms.row(i).swap(_residualNorms.row(iMin + i));
    }
  }

  /**
   * This really is just printing the eigenvalues and residual norms.
   * Any normalization cannot be done without recalculating all
   * doubles amplitudes, so that is handled by the CC2 framework.
   */
  if (_method == Options::LR_METHOD::ADC2) {
    printf("\n       ADC(2) excitation energies      \n");
  }
  else if (_method == Options::LR_METHOD::CISDINF) {
    printf("\n      CIS(Dinf) excitation energies    \n");
  }
  else if (_method == Options::LR_METHOD::CISD) {
    printf("\n        CIS(D) energy corrections      \n");
  }
  else {
    printf("\n         CC2 excitation energies       \n");
  }
  printf(" ------------------------------------------\n");
  printf("    root        eigenvalue       res norm  \n");
  printf("            (a.u.)      (eV)               \n");
  printf(" ------------------------------------------\n");
  for (unsigned i = 0; i < _nEigen; ++i) {
    printf("%6i %12.7f %9.5f %11.3e\n", i + 1, _eigenvalues(i), _eigenvalues(i) * HARTREE_TO_EV, _residualNorms(i));
  }
  printf("\n");
} /* this->postProcessing() */

} /* namespace Serenity */
