/**
 * @file ResponseSolver.cpp
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
#include "postHF/LRSCF/Tools/ResponseSolver.h"
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "parameters/Constants.h"
#include "settings/LRSCFOptions.h"

namespace Serenity {

ResponseSolver::ResponseSolver(const Eigen::VectorXd& diagonal, double convergenceCriterion, unsigned maxIterations,
                               unsigned maxSubspaceDimension, std::vector<double>& frequencies, double damping,
                               std::vector<Eigen::MatrixXd> ppmq, SigmaCalculator sigmaCalculator,
                               std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess,
                               std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk)
  : IterativeSolver(ppmq.size(), diagonal.size(), 3 * frequencies.size(), diagonal, convergenceCriterion, maxIterations,
                    maxSubspaceDimension, sigmaCalculator, writeToDisk, initialGuess),
    _frequencies(frequencies),
    _damping(damping),
    _damped(_damping == 0.0 ? false : true),
    _nFreqs(_frequencies.size()),
    _ppmq(ppmq) {
  this->initialize();
}

void ResponseSolver::initialize() {
  printHeader("Response Solver");
  if (_damped) {
    printf("%26s %16s %6.4f\n", "Damping parameter (eV)", ":", _damping * HARTREE_TO_EV);
  }

  // Set dimensions
  _expansionVectors.resize(_nSets);
  _guessVectors.resize(_nSets);
  _sigmaVectors.resize(_nSets);
  _eigenvectors.resize(_nSets);
  _residualVectors.resize(_nSets);
  _correctionVectors.resize(_nSets);
  _residualNorms = Eigen::VectorXd::Zero(_nEigen);
  _freqConverged = Eigen::VectorXi::Zero(_nFreqs);

  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _guessVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    _eigenvectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    _residualVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    _correctionVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
  }

  _eigenvalues.resize(_nFreqs * 3);
  for (unsigned iFreq = 0; iFreq < _nFreqs * 3; ++iFreq) {
    _eigenvalues(iFreq) = _frequencies[iFreq / 3];
  }

  // Initialize guessvectors
  this->seed();

  // Orthonormalize guessvectors and get rid of linear dependencies.
  Orthogonalization::modifiedGramSchmidtLinDep(_guessVectors, 1e-3);
} /* this->initialize() */

void ResponseSolver::iterate() {
  ++_nIter;

  if (_nIter == 1) {
    _itStart = std::chrono::steady_clock::now();
    _sigmaVectors = (*_sigmaCalculator(_guessVectors));
  }

  // Setup reduced system of linear equations problem
  _subDim = _guessVectors[0].cols();
  _metric = Eigen::MatrixXd::Zero(_nSets * _subDim, _nSets * _subDim);
  _subspaceMatrix = Eigen::MatrixXd::Zero(_nSets * _subDim, _nSets * _subDim);
  _subspaceEigenvectors.resize(_nSets * _subDim, _nEigen);
  Eigen::MatrixXd rhs(_nSets * _subDim, 3);

  // Construct frequency independent part of subspace matrix
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _expansionVectors[iSet].resize(_subDim, 3 * _nFreqs);
    _subspaceMatrix.block(iSet * _subDim, iSet * _subDim, _subDim, _subDim) =
        _guessVectors[iSet].transpose() * _sigmaVectors[iSet];
    _metric.block(iSet * _subDim, iSet * _subDim, _subDim, _subDim) = _guessVectors[iSet].transpose() * _guessVectors[iSet];
    rhs.block(iSet * _subDim, 0, _subDim, 3) = _guessVectors[iSet].transpose() * _ppmq[iSet];
    if (_damped) {
      _subspaceMatrix.block(iSet * _subDim, (3 - iSet) * _subDim, _subDim, _subDim) =
          (iSet < 2 ? 1. : -1.) * _damping * _guessVectors[iSet].transpose() * _guessVectors[3 - iSet];
    }
  }
  Eigen::MatrixXd cond = _metric.diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
  rhs = (cond * rhs).eval();

  Timings::takeTime("LRSCF -  RSolver: LES Solving");

  // Parallelize over frequencies
#pragma omp parallel for schedule(dynamic)
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    // A matrix for each thread
    Eigen::MatrixXd lhs = _subspaceMatrix;

    if (_freqConverged(iFreq)) {
      continue;
    }
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      unsigned iSetc = (_nSets == 1) ? 0 : (iSet / 2 * 2 + (iSet + 1) % 2);
      lhs.block(iSet * _subDim, iSetc * _subDim, _subDim, _subDim) -=
          _frequencies[iFreq] * _guessVectors[iSet].transpose() * _guessVectors[iSetc];
    } /* Loop over sets */

    // Solve subspace linear system of equations
    _subspaceEigenvectors.middleCols(3 * iFreq, 3) = cond * (cond * lhs * cond).householderQr().solve(rhs);

    // Ritz vectors
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      _expansionVectors[iSet].middleCols(3 * iFreq, 3) = _subspaceEigenvectors.block(iSet * _subDim, 3 * iFreq, _subDim, 3);
      _eigenvectors[iSet].middleCols(3 * iFreq, 3) = _guessVectors[iSet] * _expansionVectors[iSet].middleCols(3 * iFreq, 3);
    }

    // Residual vectors
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      unsigned iSetc = (_nSets == 1) ? 0 : (iSet / 2 * 2 + (iSet + 1) % 2);
      _residualVectors[iSet].middleCols(3 * iFreq, 3) =
          _sigmaVectors[iSet] * _expansionVectors[iSet].middleCols(3 * iFreq, 3);
      _residualVectors[iSet].middleCols(3 * iFreq, 3) -=
          _ppmq[iSet] + _frequencies[iFreq] * _eigenvectors[iSetc].middleCols(3 * iFreq, 3);
      if (_damped) {
        _residualVectors[iSet].middleCols(3 * iFreq, 3) +=
            (iSet < 2 ? 1. : -1.) * _damping * _eigenvectors[3 - iSet].middleCols(3 * iFreq, 3);
      }
    } /* Loop over sets */
  }   /* Loop over frequencies */
  Timings::timeTaken("LRSCF -  RSolver: LES Solving");

  // Residualnorms for each root
  _residualNorms.setZero();
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _residualNorms += _residualVectors[iSet].colwise().norm() / std::sqrt(_nSets);
  }

  // Perform convergence check frequency wise
  _nConverged = 0;
  for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
    unsigned counter = 0;
    for (unsigned iDim = 0; iDim < 3; ++iDim) {
      if (_residualNorms(3 * iFreq + iDim) < _convergenceCriterion) {
        ++_nConverged;
        ++counter;
      }
    }
    if (counter == 3) {
      _freqConverged(iFreq) = 1;
    }
  }

  _itEnd = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();

  unsigned iMax;
  _residualNorms.maxCoeff(&iMax);
  // Print Iteration Data
  if (_nIter == 1) {
    printTableHead("\n   it.    dimension    time (min)    converged     max. norm");
  }
  printf("%5i %10i %14.3f %11i %13.2e  (%2i)\n", _nIter, (int)_subDim, duration / 60.0, _nConverged,
         _residualNorms.maxCoeff(), iMax + 1);

  _itStart = std::chrono::steady_clock::now();

  // Skip the rest of this iteration if all solution vectors converged.
  if (_nConverged == _nEigen || _nIter == _maxIterations) {
    // This is done for the solver to be able to be conveniently restarted
    // from the last solution, for example, if switching to a larger grid.
    _guessVectors = _eigenvectors;

    // This is done so that in the case of switching to a larger grid, this solver
    // doesn't just skip each frequency because it thinks it is converged.
    _freqConverged = Eigen::VectorXi::Zero(_nFreqs);

    return;
  }

  // Correction vectors
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _correctionVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen - _nConverged);
  }

  // Precondition correction vectors
  this->precondition();

  // Expand subspace
  this->expandSubspace();
} /* this->iterate() */

void ResponseSolver::seed() {
  // Precondition guessvectors using orbital energy differences and the given frequency
  for (unsigned i = 0; i < _nEigen; ++i) {
    // Get corresponding frequency (it is more stable to use a finite value for omega in the static limit)
    double frequency = (_frequencies[i / 3] == 0) ? 1e-9 : _frequencies[i / 3];
    for (unsigned ia = 0; ia < _nDimension; ++ia) {
      // Invert upper left block
      double a_j = 1. / (_diagonal(ia) - frequency / _diagonal(ia) * frequency);
      double b_j = a_j * frequency / _diagonal(ia);

      if (_nSets == 1) {
        _guessVectors[0](ia, i) = 1 / (frequency - _diagonal(ia)) * _ppmq[0](ia, i % 3);
      }
      else if (_nSets == 2) {
        _guessVectors[0](ia, i) = a_j * _ppmq[0](ia, i % 3) + b_j * _ppmq[1](ia, i % 3);
        _guessVectors[1](ia, i) = b_j * _ppmq[0](ia, i % 3) + a_j * _ppmq[1](ia, i % 3);
      }
      else if (_nSets == 4) {
        double c_j = _diagonal(ia) + a_j * _damping * _damping;
        double d_j = -frequency + b_j * _damping * _damping;

        double e_j = 1. / (c_j - d_j / c_j * d_j);
        double f_j = -e_j * d_j / c_j;

        double g_j = a_j * f_j * _damping + b_j * e_j * _damping;
        double h_j = a_j * e_j * _damping + b_j * f_j * _damping;

        _guessVectors[0](ia, i) =
            e_j * _ppmq[0](ia, i % 3) + f_j * _ppmq[1](ia, i % 3) - g_j * _ppmq[2](ia, i % 3) - h_j * _ppmq[3](ia, i % 3);
        _guessVectors[1](ia, i) =
            f_j * _ppmq[0](ia, i % 3) + e_j * _ppmq[1](ia, i % 3) - h_j * _ppmq[2](ia, i % 3) - g_j * _ppmq[3](ia, i % 3);
        _guessVectors[2](ia, i) =
            g_j * _ppmq[0](ia, i % 3) + h_j * _ppmq[1](ia, i % 3) + e_j * _ppmq[2](ia, i % 3) + f_j * _ppmq[3](ia, i % 3);
        _guessVectors[3](ia, i) =
            h_j * _ppmq[0](ia, i % 3) + g_j * _ppmq[1](ia, i % 3) + f_j * _ppmq[2](ia, i % 3) + e_j * _ppmq[3](ia, i % 3);
      }
    }
  } /* End of initial guessVectors */
} /* this->seed() */

void ResponseSolver::precondition() {
  // Precondition guessvectors using orbital energy differences and the given frequency
  int index = -1;
  for (unsigned i = 0; i < _nEigen; ++i) {
    // Get corresponding frequency (it is more stable to use a finite value for omega in the static limit)
    double frequency = (_frequencies[i / 3] == 0) ? 1e-9 : _frequencies[i / 3];
    if (_residualNorms(i) < _convergenceCriterion) {
      continue;
    }
    index += 1;
    for (unsigned ia = 0; ia < _nDimension; ++ia) {
      // Invert upper left block
      double a_j = 1. / (_diagonal(ia) - frequency / _diagonal(ia) * frequency);
      double b_j = a_j * frequency / _diagonal(ia);

      if (_nSets == 1) {
        _correctionVectors[0](ia, index) = 1 / (frequency - _diagonal(ia)) * _residualVectors[0](ia, i);
      }
      else if (_nSets == 2) {
        _correctionVectors[0](ia, index) = a_j * _residualVectors[0](ia, i) + b_j * _residualVectors[1](ia, i);
        _correctionVectors[1](ia, index) = b_j * _residualVectors[0](ia, i) + a_j * _residualVectors[1](ia, i);
      }
      else if (_nSets == 4) {
        double c_j = _diagonal(ia) + a_j * _damping * _damping;
        double d_j = -frequency + b_j * _damping * _damping;

        double e_j = 1. / (c_j - d_j / c_j * d_j);
        double f_j = -e_j * d_j / c_j;

        double g_j = a_j * f_j * _damping + b_j * e_j * _damping;
        double h_j = a_j * e_j * _damping + b_j * f_j * _damping;

        _correctionVectors[0](ia, index) = e_j * _residualVectors[0](ia, i) + f_j * _residualVectors[1](ia, i) -
                                           g_j * _residualVectors[2](ia, i) - h_j * _residualVectors[3](ia, i);
        _correctionVectors[1](ia, index) = f_j * _residualVectors[0](ia, i) + e_j * _residualVectors[1](ia, i) -
                                           h_j * _residualVectors[2](ia, i) - g_j * _residualVectors[3](ia, i);
        _correctionVectors[2](ia, index) = g_j * _residualVectors[0](ia, i) + h_j * _residualVectors[1](ia, i) +
                                           e_j * _residualVectors[2](ia, i) + f_j * _residualVectors[3](ia, i);
        _correctionVectors[3](ia, index) = h_j * _residualVectors[0](ia, i) + g_j * _residualVectors[1](ia, i) +
                                           f_j * _residualVectors[2](ia, i) + e_j * _residualVectors[3](ia, i);
      }
    }
  }
} /* this->precondition */

} /* namespace Serenity */
