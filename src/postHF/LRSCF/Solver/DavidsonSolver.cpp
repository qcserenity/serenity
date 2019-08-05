/**
 * @file DavidsonSolver.cpp
 *
 * @date Aug 31, 2016
 * @author Michael Boeckers
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
#include "postHF/LRSCF/Solver/DavidsonSolver.h"
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <iomanip>


namespace Serenity {

DavidsonSolver::DavidsonSolver(
    unsigned int nDimension,
    unsigned int nEigen,
    Eigen::VectorXd& diagonalElements,
    double convergenceCriterion,
    unsigned int nMaxIterations,
    std::function<Eigen::MatrixXd (
        Eigen::MatrixXd& guessVectors)> sigmaVectorCalculator,
    unsigned int nMaxSubspaceDimension,
    std::string loadPath,
    std::string filePath,
    std::string id):
      IterativeEigenvalueSolver(
          nDimension,
          nEigen,
          diagonalElements,
          convergenceCriterion,
          nMaxIterations,
          loadPath,
          filePath,
          id),
      _sigmaVectorCalculator(sigmaVectorCalculator),
      _nMaxSubspaceDimension(nMaxSubspaceDimension){
  initialize();
}

void DavidsonSolver::iterate() {
  //set up subspace matrix
  _subspaceMatrix = _guessVectors[0].transpose() * _sigmaVectors[0];

  //Compute eigenpairs of reduced eigenvalue problem.
  //Use qr eigenvalue solver of eigen library.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(_guessVectors[0].cols());
  eigenvalueSolver.compute(_subspaceMatrix);
  _subspaceEigenvalues = eigenvalueSolver.eigenvalues();
  _subspaceEigenvectors = eigenvalueSolver.eigenvectors();
  _eigenvalues = _subspaceEigenvalues.block(0,0,_nEigen,1);

  //calculate current approximation to eigenvectors
  computeTrialVectors(
      _eigenvectors[0],
      _subspaceEigenvectors,
      _guessVectors[0]);

  //calculate residuals
  computeResidualVectors(
      _residualVectors[0],
      _sigmaVectors[0],
      _guessVectors[0],
      _subspaceEigenvectors,
      _subspaceEigenvalues);
  checkForConvergence(_residualVectors);

  //calculate correction vectors
  Eigen::MatrixXd correctionVectors;
  computeCorrectionVectors(
      correctionVectors,
      _residualVectors[0]);

  //update test vectors
  updateGuessVectors(
      correctionVectors,
      _guessVectors[0]);

  //perform subspace collapse if dimension above threshold
  subspaceCollapse(
      correctionVectors,
      _guessVectors[0],
      _sigmaVectors[0],
      _eigenvectors[0]);

  //compute new sigma vectors
  Eigen::MatrixXd newSigmaVectors;
  newSigmaVectors = _sigmaVectorCalculator(correctionVectors);

  //update sigma vectors
  unsigned int oldDimension = _sigmaVectors[0].cols();
  unsigned int _nEigen = oldDimension + newSigmaVectors.cols();
  _sigmaVectors[0].conservativeResize(_nDimension,_nEigen);
  _sigmaVectors[0].block(0,oldDimension,_nDimension,newSigmaVectors.cols()) = newSigmaVectors;

  //print iteration data
  printIterationData(_residualVectors,_eigenvalues);
}

void DavidsonSolver::postProcessing() {
  return;
}


void DavidsonSolver::computeCorrectionVectors(
    Eigen::MatrixXd& correctionVectors,
    Eigen::MatrixXd& residualVectors) {
  correctionVectors.resize(_nDimension, _nEigen - _nConverged);
  int iCount = -1;
  for (int i = 0; i < residualVectors.cols(); ++i) {
    if (residualVectors.col(i).norm() < _convergenceCriterion) {
      continue;
    }
    iCount += 1;
    Eigen::VectorXd diagonalPreconditioner = _subspaceEigenvalues(i) *
        Eigen::VectorXd::Ones(_nDimension) - _diagonalElements;
    for (unsigned int j = 0; j < _nDimension; ++j) {
      correctionVectors(j, iCount) = residualVectors(j, i) / diagonalPreconditioner(j);
    }
  }
  correctionVectors.colwise().normalize();
}

void DavidsonSolver::subspaceCollapse(
    Eigen::MatrixXd& correctionVectors,
    Eigen::MatrixXd& guessVectors,
    Eigen::MatrixXd& sigmaVectors,
    Eigen::MatrixXd& eigenvectors) {
  if (guessVectors.cols() >= _nMaxSubspaceDimension) {
    std::printf("perform subspace collapse\n");
    guessVectors= eigenvectors;
    //Must calculate new sigma vectors for this new set of guess vectors.
    //Since we do not want to neglect the current correction vectors, we
    //also append them.
    guessVectors.conservativeResize(_nDimension,_nEigen + correctionVectors.cols());
    guessVectors.block(0,_nEigen,_nDimension,correctionVectors.cols()) = correctionVectors;
    Eigen::MatrixXd newSigmas;
    computeTrialVectors(newSigmas,_subspaceEigenvectors,sigmaVectors);
    sigmaVectors = newSigmas;
  }
}

void DavidsonSolver::updateGuessVectors(
    Eigen::MatrixXd& correctionVectors,
    Eigen::MatrixXd& guessVectors) {
  int oldSubspaceDimension = guessVectors.cols();
  //orthogonalize correction vectors agains set of guess vectors
  Eigen::MatrixXd testMatrix(_nDimension,oldSubspaceDimension + correctionVectors.cols());
  testMatrix.block(0,0,guessVectors.rows(),oldSubspaceDimension) = guessVectors;
  testMatrix.block(0,oldSubspaceDimension,guessVectors.rows(),correctionVectors.cols()) = correctionVectors;
  bool orthogonal = false;
  unsigned int nIter = 0;
  //iterative orthogonalization
  do {
    nIter += 1;
    Orthogonalization::modifiedGramSchmidt (testMatrix);
    //test orthogonality
    Eigen::MatrixXd thisShouldBeIdentity = guessVectors.transpose() * guessVectors;
    double offDiagonalSum = 0;
    for (int i = 0; i < thisShouldBeIdentity.rows(); ++i) {
      for (int j = 0; j < thisShouldBeIdentity.cols(); ++j) {
        if (i != j) {
          offDiagonalSum += std::abs(thisShouldBeIdentity(i,j));
        }
      }
    }
    if (offDiagonalSum < 1.0e-8) {
      orthogonal = true;
//      print((std::string) "\n Orthogonality ensured after " + nIter + " iterations" );
    }
    if ((nIter > 10 and offDiagonalSum > 1.0e-3) or
        (nIter >= 100)) {
      throw SerenityError("could not orthonormalize test vectors.");
    }
  } while (orthogonal == false);
  //remove correction vectors with Schmidt-norm lower than threshold
  int nAppend = 0;
  for (int i = 0; i < correctionVectors.cols(); ++i) {
    if (testMatrix.col(oldSubspaceDimension + i).norm() >= 1.0e-4) {
      nAppend += 1;
      correctionVectors.col(nAppend -1) = testMatrix.col(oldSubspaceDimension + i);
    }
  }
  //append correction vectors and normalize
  guessVectors.resize(_nDimension,oldSubspaceDimension + nAppend);
  guessVectors.block(0,0,_nDimension,oldSubspaceDimension) = testMatrix.block(0,0,_nDimension,oldSubspaceDimension);
  guessVectors.block(0,oldSubspaceDimension,_nDimension,nAppend) = correctionVectors.block(0,0,_nDimension,nAppend);
  guessVectors.colwise().normalize();
  //modify correction vectors
  correctionVectors.resize(_nDimension,nAppend);
  correctionVectors = guessVectors.block(0,oldSubspaceDimension,_nDimension,nAppend);
}


void DavidsonSolver::initialize() {
  //print header
  printHeader("Davidson eigenvalue solver");
  print((std::string)"Maximum subspace dimension          : " + _nMaxSubspaceDimension);
  //Only one set of guess vectors is used
  _guessVectors.resize(1);
  _residualVectors.resize(1);
  _sigmaVectors.resize(1);
  _eigenvectors.resize(1);
  if (!initializeFromHDF5()) {
    //calculate sigma vectors
    computeInitialGuessVectors(_guessVectors[0],_diagonalElements);
    _sigmaVectors[0] = _sigmaVectorCalculator(_guessVectors[0]);
  } else {
    std::cout << "RESTART \n";
    if (_guessVectors[0].cols() < _nEigen) {
      //ToDo: Reuse guess vectors from HDF5 file...
      computeInitialGuessVectors(_guessVectors[0],_diagonalElements);
      _sigmaVectors[0] = _sigmaVectorCalculator(_guessVectors[0]);
    }
  }
}



} /* namespace Serenity */
