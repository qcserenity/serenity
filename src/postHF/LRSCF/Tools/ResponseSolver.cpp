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

ResponseSolver::ResponseSolver(
    unsigned nDimension, unsigned nEigen, Eigen::VectorXd& diagonal, double convergenceCriterion, unsigned maxIterations,
    unsigned maxSubspaceDimension, std::vector<double>& frequencies, double damping, Options::GAUGE gauge,
    const Eigen::MatrixXd& lengths, const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& magnetics,
    std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)> sigmaCalculator,
    std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess)
  : IterativeSolver(damping == 0.0 ? 2 : 4, nDimension, 0 * nEigen + 3 * frequencies.size(), diagonal,
                    convergenceCriterion, maxIterations, maxSubspaceDimension, sigmaCalculator, initialGuess),
    _frequencies(frequencies),
    _damping(damping),
    _gauge(gauge),
    _lengths(lengths),
    _velocities(velocities),
    _magnetics(magnetics) {
  this->initialize();
}

void ResponseSolver::initialize() {
  printHeader("Response Solver");
  if (_damped)
    printf("%26s %16s %6.4f\n", "Damping Parameter (eV)", ":", _damping * HARTREE_TO_EV);

  // Set class member variables
  _damped = _damping == 0.0 ? false : true;
  _nFreqs = _frequencies.size();

  // Set dimensions
  _ppmq.resize(_nSets);
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
    _ppmq[iSet] = Eigen::MatrixXd::Zero(_nDimension, 3);
  }

  /**
   * The right-hand side of the linear equation system. Correctly handling the signs in a
   * systematic fashion (and not just correcting them afterwards) led to a lot of confusion
   * on my side and it is NOWHERE explicitly explained in the literature, so I will
   * do that here. The general TDDFT response equation reads
   *
   *                         (A B)(X)-w(Y) = -Q
   *                         (B A)(Y)+w(Y) = -R
   *
   * so the the outer minus (from below) comes from here. Elements in P and Q are defined such that:
   *
   *                         Q_ia = <i|V|a>,
   *                         R_ia = P_ia^* = <a|V|i>
   *
   * and V describes the perturbation. We want to obtain the polarizability as the reaction
   * of the electric dipole to an electric-dipole perturbation. In response theory, we can
   * describe this with the linear-response function <<A;V>>, where A describes the observable
   * (A = mu) and V describes the perturbing field. Here we note that we must substitute
   * V = -mu F (F is the merely the amplitude of the field). We thus insert the operator -mu
   * in the response function and the above property integrals (the inner minus comes from here).
   *
   * We solve for -(Q+R) (stored in _ppmq[0]) and -(Q-R) (stored in _ppmq[1]) here. We thus must note that
   * mu (length) is real and p (velocity) is imaginary. Since both operators are Hermitian,
   * we find in the case of V = -mu
   *
   *                         -(Q+R) = -(-2mu) and -(Q-R) = 0 in the length case (since <i|mu|a> = <a|mu|i>)
   *
   * and in the case of V = -p, we find
   *
   *                         -(Q+R) = 0 and -(Q-R) = -(-2p) in the velocity case (since <i|p|a> = -<a|p|i>).
   */
  if (_gauge == Options::GAUGE::LENGTH) {
    _ppmq[0] = -(-2 * _lengths);
  }
  else if (_gauge == Options::GAUGE::VELOCITY) {
    _ppmq[1] = -(-2 * _velocities);
  }

  // Initialize guessvectors
  this->seed();

  // Orthonormalize guessvectors and get rid of linear dependencies
  // Note: if several frequencies close to each other were selected
  //       the initial guessvectors are prone to linear dependencies
  //       due to the seeding procedure
  Orthogonalization::modifiedGramSchmidtLinDep(_guessVectors, 1e-4);

  // Calculate initial sigma vectors
  _sigmaVectors = (*_sigmaCalculator(_guessVectors));
} /* this->initialize() */

void ResponseSolver::iterate() {
  // Setup reduced system of linear equations problem
  _subDim = _guessVectors[0].cols();
  _subspaceMatrix = Eigen::MatrixXd::Zero(_nSets * _subDim, _nSets * _subDim);
  _subspaceEigenvectors.resize(_nSets * _subDim, _nEigen);
  Eigen::MatrixXd rhs(_nSets * _subDim, 3);

  // Construct frequency independent part of subspace matrix
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _expansionVectors[iSet].resize(_subDim, 3 * _nFreqs);
    _subspaceMatrix.block(iSet * _subDim, iSet * _subDim, _subDim, _subDim) =
        _guessVectors[iSet].transpose() * _sigmaVectors[iSet];
    rhs.block(iSet * _subDim, 0, _subDim, 3) = _guessVectors[iSet].transpose() * _ppmq[iSet];
    if (_damped) {
      _subspaceMatrix.block(iSet * _subDim, (3 - iSet) * _subDim, _subDim, _subDim) =
          (iSet < 2 ? 1. : -1.) * _damping * _guessVectors[iSet].transpose() * _guessVectors[3 - iSet];
    }
  }

  Timings::takeTime("LRSCF -  RSolver: LES Solving");

  // Parallelize over frequencies
  // Safety first
  Eigen::setNbThreads(1);
#pragma omp parallel
  {
    // A matrix for each thread
    Eigen::MatrixXd lhs = _subspaceMatrix;

#pragma omp for schedule(dynamic)
    for (unsigned iFreq = 0; iFreq < _nFreqs; ++iFreq) {
      if (_freqConverged(iFreq))
        continue;
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        unsigned iSetc = iSet / 2 * 2 + (iSet + 1) % 2;
        lhs.block(iSet * _subDim, iSetc * _subDim, _subDim, _subDim) =
            -_frequencies[iFreq] * _guessVectors[iSet].transpose() * _guessVectors[iSetc];
      } /* Loop over sets */

      // Solve subspace linear system of equations
      _subspaceEigenvectors.middleCols(3 * iFreq, 3) = lhs.householderQr().solve(rhs);

      // Ritz vectors
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        _expansionVectors[iSet].middleCols(3 * iFreq, 3) = _subspaceEigenvectors.block(iSet * _subDim, 3 * iFreq, _subDim, 3);
        _eigenvectors[iSet].middleCols(3 * iFreq, 3) = _guessVectors[iSet] * _expansionVectors[iSet].middleCols(3 * iFreq, 3);
      }

      // Residual vectors
      for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
        unsigned iSetc = iSet / 2 * 2 + (iSet + 1) % 2;
        _residualVectors[iSet].middleCols(3 * iFreq, 3) =
            _sigmaVectors[iSet] * _expansionVectors[iSet].middleCols(3 * iFreq, 3) - _ppmq[iSet];
        _residualVectors[iSet].middleCols(3 * iFreq, 3) -=
            _frequencies[iFreq] * _eigenvectors[iSetc].middleCols(3 * iFreq, 3);
        if (_damped) {
          _residualVectors[iSet].middleCols(3 * iFreq, 3) +=
              (iSet < 2 ? 1. : -1.) * _damping * _eigenvectors[3 - iSet].middleCols(3 * iFreq, 3);
        }
      } /* Loop over sets */
    }   /* Loop over frequencies */
  }     /* pragma omp for */
  Eigen::setNbThreads(0);
  Timings::timeTaken("LRSCF -  RSolver: LES Solving");

  // Residualnorms for each root
  _residualNorms.setZero();
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _residualNorms += _residualVectors[iSet].colwise().norm() / std::sqrt(_nSets);
  }

  // Convergence check
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
    if (counter == 3)
      _freqConverged(iFreq) = 1;
  }

  // Print Iteration Data
  if (_nIter == 1)
    printTableHead("\n   It.    Dimension    Converged    Max. Norm");
  printf("%5i %11i %12i %14.2e\n", _nIter, (int)_subDim, _nConverged, _residualNorms.maxCoeff());

  // Skip the rest of this iteration if all solution vectors converged
  if (_nConverged == _nEigen)
    return;

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
    for (unsigned j = 0; j < _nDimension; ++j) {
      // Invert upper left block
      double a_j = 1. / (_diagonal(j) - frequency / _diagonal(j) * frequency);
      double b_j = a_j * frequency / _diagonal(j);

      if (!_damped) {
        _guessVectors[0](j, i) = a_j * _ppmq[0](j, i % 3) + b_j * _ppmq[1](j, i % 3);
        _guessVectors[1](j, i) = b_j * _ppmq[0](j, i % 3) + a_j * _ppmq[1](j, i % 3);
      }
      else if (_damped) {
        double c_j = _diagonal(j) + a_j * _damping * _damping;
        double d_j = -frequency + b_j * _damping * _damping;

        double e_j = 1. / (c_j - d_j / c_j * d_j);
        double f_j = -e_j * d_j / c_j;

        double g_j = a_j * f_j * _damping + b_j * e_j * _damping;
        double h_j = a_j * e_j * _damping + b_j * f_j * _damping;

        _guessVectors[0](j, i) =
            e_j * _ppmq[0](j, i % 3) + f_j * _ppmq[1](j, i % 3) - g_j * _ppmq[2](j, i % 3) - h_j * _ppmq[3](j, i % 3);
        _guessVectors[1](j, i) =
            f_j * _ppmq[0](j, i % 3) + e_j * _ppmq[1](j, i % 3) - h_j * _ppmq[2](j, i % 3) - g_j * _ppmq[3](j, i % 3);
        _guessVectors[2](j, i) =
            g_j * _ppmq[0](j, i % 3) + h_j * _ppmq[1](j, i % 3) + e_j * _ppmq[2](j, i % 3) + f_j * _ppmq[3](j, i % 3);
        _guessVectors[3](j, i) =
            h_j * _ppmq[0](j, i % 3) + g_j * _ppmq[1](j, i % 3) + f_j * _ppmq[2](j, i % 3) + e_j * _ppmq[3](j, i % 3);
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
    if (_residualNorms(i) < _convergenceCriterion)
      continue;
    index += 1;
    for (unsigned j = 0; j < _nDimension; ++j) {
      // Invert upper left block
      double a_j = 1. / (_diagonal(j) - frequency / _diagonal(j) * frequency);
      double b_j = a_j * frequency / _diagonal(j);

      if (!_damped) {
        _correctionVectors[0](j, index) = a_j * _residualVectors[0](j, i) + b_j * _residualVectors[1](j, i);
        _correctionVectors[1](j, index) = b_j * _residualVectors[0](j, i) + a_j * _residualVectors[1](j, i);
      }
      else if (_damped) {
        double c_j = _diagonal(j) + a_j * _damping * _damping;
        double d_j = -frequency + b_j * _damping * _damping;

        double e_j = 1. / (c_j - d_j / c_j * d_j);
        double f_j = -e_j * d_j / c_j;

        double g_j = a_j * f_j * _damping + b_j * e_j * _damping;
        double h_j = a_j * e_j * _damping + b_j * f_j * _damping;

        _correctionVectors[0](j, index) = e_j * _residualVectors[0](j, i) + f_j * _residualVectors[1](j, i) -
                                          g_j * _residualVectors[2](j, i) - h_j * _residualVectors[3](j, i);
        _correctionVectors[1](j, index) = f_j * _residualVectors[0](j, i) + e_j * _residualVectors[1](j, i) -
                                          h_j * _residualVectors[2](j, i) - g_j * _residualVectors[3](j, i);
        _correctionVectors[2](j, index) = g_j * _residualVectors[0](j, i) + h_j * _residualVectors[1](j, i) +
                                          e_j * _residualVectors[2](j, i) + f_j * _residualVectors[3](j, i);
        _correctionVectors[3](j, index) = h_j * _residualVectors[0](j, i) + g_j * _residualVectors[1](j, i) +
                                          f_j * _residualVectors[2](j, i) + e_j * _residualVectors[3](j, i);
      }
    }
  }
} /* this->precondition */
} /* namespace Serenity */
