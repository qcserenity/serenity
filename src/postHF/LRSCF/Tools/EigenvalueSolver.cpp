/**
 * @file EigenvalueSolver.cpp
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
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/SerenityError.h"
#include "parameters/Constants.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <fstream>
#include <iomanip>

namespace Serenity {

EigenvalueSolver::EigenvalueSolver(
    unsigned nDimension, unsigned nEigen, const Eigen::VectorXd& diagonal, double convergenceCriterion,
    unsigned maxIterations, unsigned maxSubspaceDimension, unsigned initialSubspace, Options::LR_METHOD algorithm,
    std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)> sigmaCalculator,
    std::shared_ptr<std::vector<Eigen::MatrixXd>> initialGuess,
    std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk)
  : IterativeSolver(algorithm == Options::LR_METHOD::TDDFT ? 2 : 1, nDimension, nEigen, diagonal, convergenceCriterion,
                    maxIterations, maxSubspaceDimension, sigmaCalculator, writeToDisk, initialGuess),
    _initialSubspace(initialSubspace),
    _algorithm(algorithm) {
  this->initialize();
}

void EigenvalueSolver::initialize() {
  printHeader("Eigenvalue Solver");
  if (_algorithm == Options::LR_METHOD::TDA) {
    printf("\n    TDA algorithm.\n");
  }
  else if (_algorithm == Options::LR_METHOD::TDDFT) {
    printf("\n    TDDFT algorithm.\n");
  }
  else {
    throw SerenityError("Cannot use linear eigenvalue solver for CC2 methods.");
  }

  // set dimensions
  _sigmaVectors.resize(_nSets);
  _guessVectors.resize(_nSets);
  _eigenvectors.resize(_nSets);
  _residualVectors.resize(_nSets);
  _expansionVectors.resize(_nSets);
  _correctionVectors.resize(_nSets);
  _residualNorms = Eigen::VectorXd::Zero(_nEigen);
  _eigenvalues = Eigen::VectorXd::Zero(_nEigen);

  // initialize all matrices
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _guessVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _initialSubspace);
    _eigenvectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    _residualVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    _correctionVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
  }

  // fill up guess vectors with unit vectors with ones
  // at the positions with the lowest orb-energy differences
  std::vector<std::pair<double, unsigned>> diagPair;
  for (unsigned ia = 0; ia < _nDimension; ++ia) {
    diagPair.push_back(std::pair<double, unsigned>(_diagonal(ia), ia));
  }
  std::stable_sort(diagPair.begin(), diagPair.end());
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    for (unsigned i = 0; i < _initialSubspace; ++i) {
      _guessVectors[iSet](diagPair[i].second, i) = 1.0;
    }
  }

  if (_initialGuess) {
    unsigned inputRoots = (*_initialGuess)[0].cols();
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      _guessVectors[iSet].rightCols(inputRoots) = (*_initialGuess)[iSet];
    }
    Orthogonalization::modifiedGramSchmidtLinDep(_guessVectors, (inputRoots < _nEigen) ? 1e-3 : 0.0);
  }
} /* this->initialize() */

void EigenvalueSolver::iterate() {
  ++_nIter;

  if (_nIter == 1) {
    _itStart = std::chrono::steady_clock::now();
    _sigmaVectors = (*_sigmaCalculator(_guessVectors));
  }

  _subDim = _guessVectors[0].cols();

  _metric = Eigen::MatrixXd::Zero(_nSets * _subDim, _nSets * _subDim);
  _subspaceMatrix = Eigen::MatrixXd::Zero(_nSets * _subDim, _nSets * _subDim);
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _subspaceMatrix.block(iSet * _subDim, iSet * _subDim, _subDim, _subDim) =
        _guessVectors[iSet].transpose() * _sigmaVectors[iSet];
    _metric.block(iSet * _subDim, iSet * _subDim, _subDim, _subDim) = _guessVectors[iSet].transpose() * _guessVectors[iSet];
  }
  Eigen::MatrixXd cond = _metric.diagonal().cwiseSqrt().cwiseInverse().asDiagonal();

  _subspaceMatrix = (cond * _subspaceMatrix * cond).eval();

  if (_nSets == 1) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(_subspaceMatrix);
    _subspaceEigenvectors = cond * subspaceSolver.eigenvectors();
    _subspaceEigenvalues = subspaceSolver.eigenvalues();

    if (_maxIterations == 1) {
      printf("\n");
      printSmallCaption("Coupling Matrix / eV");
      for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
        printf(" ");
        for (unsigned jEigen = 0; jEigen < _nEigen; ++jEigen) {
          printf("%12.3e", _subspaceMatrix(iEigen, jEigen) * HARTREE_TO_EV);
        }
        printf("\n");
      }
      printf("\n");
    }
  }
  else {
    _metric.setZero();
    _metric.topRightCorner(_subDim, _subDim) = _guessVectors[0].transpose() * _guessVectors[1];
    _metric.bottomLeftCorner(_subDim, _subDim) = _guessVectors[1].transpose() * _guessVectors[0];
    _metric = (cond * _metric * cond).eval();

    Eigen::FullPivLU<Eigen::MatrixXd> lu(_metric);
    if (lu.isInvertible()) {
      _subspaceMatrix = (lu.inverse() * _subspaceMatrix).eval();
    }
    else {
      OutputControl::dOut << "  Subspace metric:\n" << _subspaceMatrix << std::endl;
      throw SerenityError("Subspace metric cannot be inverted.");
    }
    Eigen::EigenSolver<Eigen::MatrixXd> subspaceSolver(_subspaceMatrix);
    _subspaceEigenvectors = cond * subspaceSolver.eigenvectors().real();
    _subspaceEigenvalues = subspaceSolver.eigenvalues().real();

    // sort in ascending order
    unsigned iMin;
    for (unsigned i = 0; i < _nSets * _subDim; ++i) {
      _subspaceEigenvalues.tail(_nSets * _subDim - i).minCoeff(&iMin);
      _subspaceEigenvalues.row(i).swap(_subspaceEigenvalues.row(iMin + i));
      _subspaceEigenvectors.col(i).swap(_subspaceEigenvectors.col(iMin + i));
    }

    _subspaceEigenvalues = _subspaceEigenvalues.tail(_subDim).eval();
    _subspaceEigenvectors = _subspaceEigenvectors.rightCols(_subDim).eval();
  }

  // if there was an initial guess, the eigenvalue solver will try to follow _nEigen roots
  // that were given to it based on the overlap of the ritz vectors with the initial guess
  if (_initialGuess) {
    Eigen::MatrixXd& rootsToFollow = (_nIter == 1) ? (*_initialGuess)[0] : _eigenvectors[0];
    Eigen::MatrixXd ritzVectors = _guessVectors[0] * _subspaceEigenvectors.topRows(_subDim);
    Eigen::MatrixXd overlap = ritzVectors.transpose() * rootsToFollow;
    // sort nRootsToFollow in descending order (wrt to overlap)
    unsigned iMax;
    for (unsigned i = 0; i < rootsToFollow.cols(); ++i) {
      overlap.col(i).cwiseAbs().maxCoeff(&iMax);
      // swap eigenvalues and expansion vectors
      _subspaceEigenvalues.row(i).swap(_subspaceEigenvalues.row(iMax));
      _subspaceEigenvectors.col(i).swap(_subspaceEigenvectors.col(iMax));
    }
    // sort the rest in ascending order to catch the lowest-lying missing roots
    if (rootsToFollow.cols() < _nEigen) {
      unsigned iMin;
      for (unsigned i = rootsToFollow.cols(); i < _subDim; ++i) {
        _subspaceEigenvalues.tail(_subDim - i).minCoeff(&iMin);
        _subspaceEigenvalues.row(i).swap(_subspaceEigenvalues.row(iMin + i));
        _subspaceEigenvectors.col(i).swap(_subspaceEigenvectors.col(iMin + i));
      }
    }
  }

  // only take as many as required by _nEigen
  _subspaceEigenvalues = _subspaceEigenvalues.head(_nEigen).eval();
  _subspaceEigenvectors = _subspaceEigenvectors.leftCols(_nEigen).eval();

  // sort one last time
  unsigned iMin;
  for (unsigned i = 0; i < _nEigen; ++i) {
    _subspaceEigenvalues.tail(_nEigen - i).minCoeff(&iMin);
    _subspaceEigenvalues.row(i).swap(_subspaceEigenvalues.row(iMin + i));
    _subspaceEigenvectors.col(i).swap(_subspaceEigenvectors.col(iMin + i));
  }

  // transfer to expansion vectors
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _expansionVectors[iSet].resize(_subDim, _nEigen);
  }
  for (unsigned i = 0; i < _nEigen; ++i) {
    _eigenvalues(i) = _subspaceEigenvalues(i);
    for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
      _expansionVectors[iSet].col(i) = _subspaceEigenvectors.middleRows(iSet * _subDim, _subDim).col(i);
    }
  }

  // ritz vectors
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _eigenvectors[iSet] = _guessVectors[iSet] * _expansionVectors[iSet];
  }

  // residual vectors and norms
  _residualNorms.setZero();
  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _residualVectors[iSet] = _sigmaVectors[iSet] * _expansionVectors[iSet] - _eigenvectors[_nSets == 1 ? 0
                                                                                           : iSet == 0 ? 1
                                                                                                       : 0] *
                                                                                 _eigenvalues.asDiagonal();
    _residualNorms += _residualVectors[iSet].colwise().norm() / std::sqrt(_nSets);
  }

  // check convergence
  _nConverged = 0;
  for (unsigned i = 0; i < _nEigen; ++i) {
    if (_residualNorms(i) < _convergenceCriterion)
      ++_nConverged;
  }

  _itEnd = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::duration<double>>(_itEnd - _itStart).count();

  unsigned iMax;
  _residualNorms.maxCoeff(&iMax);
  // print iteration data
  if (_nIter == 1) {
    printTableHead("\n   it.    dimension    time (min)    converged     max. norm");
  }
  printf("%5i %10i %14.3f %11i %13.2e  (%2i)\n", _nIter, (int)_subDim, duration / 60.0, _nConverged,
         _residualNorms.maxCoeff(), iMax + 1);

  _itStart = std::chrono::steady_clock::now();

  if (_nConverged == _nEigen || _nIter == _maxIterations) {
    // This is done for the solver to be able to be conveniently restarted
    // from the last solution, for example, if switching to a larger grid.
    _guessVectors = _eigenvectors;
    return;
  }

  for (unsigned iSet = 0; iSet < _nSets; ++iSet) {
    _correctionVectors[iSet].resize(_nDimension, _nEigen - _nConverged);
    _correctionVectors[iSet].setZero();
  }

  int iCount = -1;
  for (unsigned j = 0; j < _nEigen; ++j) {
    if (_residualNorms(j) < _convergenceCriterion)
      continue;
    ++iCount;
    for (unsigned ia = 0; ia < _nDimension; ++ia) {
      if (_algorithm == Options::LR_METHOD::TDA) {
        double aj = 1 / (_eigenvalues(j) - _diagonal(ia));
        _correctionVectors[0](ia, iCount) = aj * _residualVectors[0](ia, j);
      }
      else {
        double aj = 1 / (_diagonal(ia) - _eigenvalues(j) / _diagonal(ia) * _eigenvalues(j));
        double bj = aj * _eigenvalues(j) / _diagonal(ia);
        _correctionVectors[0](ia, iCount) = aj * _residualVectors[0](ia, j) + bj * _residualVectors[1](ia, j);
        _correctionVectors[1](ia, iCount) = bj * _residualVectors[0](ia, j) + aj * _residualVectors[1](ia, j);
      }
    }
  }

  this->expandSubspace();

  this->_writeToDisk(_eigenvectors, _eigenvalues);
} /* this->iterate() */

void EigenvalueSolver::postProcessing() {
  // Normalize eigenvectors.
  if (_algorithm == Options::LR_METHOD::TDDFT) {
    // Normalize <X+Y|X-Y>.
    for (unsigned iEigen = 0; iEigen < _nEigen; ++iEigen) {
      double clr = _eigenvectors[1].col(iEigen).dot(_eigenvectors[0].col(iEigen));
      _eigenvectors[0].col(iEigen) *= 1.0 / std::sqrt(clr);
      _eigenvectors[1].col(iEigen) *= 1.0 / std::sqrt(clr);
    }
  }

  // print eigenvalues with norms
  if (_algorithm == Options::LR_METHOD::TDA) {
    printf("         TDA excitation energies       \n");
  }
  else {
    printf("        TDDFT excitation energies      \n");
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
