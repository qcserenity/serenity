/**
 * @file   KDavidsonEigenvalueSolver.cpp
 *
 * @date   Mar 4, 2017
 * @author M. Boeckers
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
#include "postHF/LRSCF/Solver/KDavidsonEigenvalueSolver.h"
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/Timing.h"


namespace Serenity {

KDavidsonEigenvalueSolver::KDavidsonEigenvalueSolver(
    unsigned int nDimension,
    unsigned int nEigen,
    Eigen::VectorXd& diagonalElements,
    double convergenceCriterion,
    unsigned int nMaxIterations,
    std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVectors)> MSigmaCalculator,
    std::function<Eigen::MatrixXd(Eigen::MatrixXd& guessVectors)> KSigmaCalculator,
    unsigned int maxSubspaceDimension,
    std::string loadPath,
    std::string filePath,
    std::string id) :
      IterativeEigenvalueSolver(
          nDimension,
          nEigen,
          diagonalElements,
          convergenceCriterion,
          nMaxIterations,
          loadPath,
          filePath,
          id),
      _MSigmaCalculator(MSigmaCalculator),
      _KSigmaCalculator(KSigmaCalculator),
      _maxSubspaceDimension(maxSubspaceDimension){
  initializeKDavidson();
  if (_maxSubspaceDimension > nDimension) _maxSubspaceDimension = _nDimension;
}


void KDavidsonEigenvalueSolver::computeCorrectionVectors(
    Eigen::MatrixXd& residualVectors,
    Eigen::MatrixXd& correctionVectors,
    VectorXd& diagonalElements) {
    correctionVectors.resize(_nDimension, _nEigen - _nConverged);
    if (correctionVectors.cols() == 0) return;
    int iCount = -1;
    for (unsigned int i = 0; i < residualVectors.cols(); ++i) {
      if (residualVectors.col(i).norm() < _convergenceCriterion) {
        continue;
      }
      iCount += 1;
//      VectorXd diagonalPreconditioner = _subspaceEigenvalues(i) * VectorXd::Ones(_nDimension) - diagonalElements;
      for (unsigned int j = 0; j < _nDimension; ++j) {
         double tmpDouble = _subspaceEigenvalues(i) - diagonalElements(j);
         if (std::abs(tmpDouble) < 1.0e-6) {
           //ensure, that correction vector does not become singular
           //ToDo: There must be more elegant solution...
           tmpDouble = 1.0e-6;
         }
         correctionVectors(j, iCount) = residualVectors(j, i) / tmpDouble;
//        correctionVectors(j, iCount) = residualVectors(j, i) / diagonalPreconditioner(j);
      }
    }
    correctionVectors.colwise().normalize();
}

void KDavidsonEigenvalueSolver::kOrthogonalize(
    Eigen::MatrixXd& S,
    Eigen::MatrixXd& KS,
    Eigen::MatrixXd& MKS) {
  //This is what actually happens and what is suggested by E. Vacharynski et. al.:
  //
  //    LLT<MatrixXd> chol(S.transpose() * KS);
  //    Eigen::MatrixXd R = chol.matrixL();
  //    Eigen::MatrixXd rinverse = R.inverse();
  //    S = S * rinverse;
  //    KS = KS * rinverse;
  //    MKS = MKS * rinverse;
  //
  //However, since the k-orthonormalization using the Cholesky decomposition is unstable when the
  //incoming matrix is ill-conditioned, we use the more stable SVQB (singular vector QB factorization).
  //Normally, this orthogonalization method requires only one or two iterations. However, in the case of
  //extremly ill-conditioned matrices S, the SVQB algorithm may use a third application of the procedure to
  //produce orthogonal vectors.
  //
  //Ref: A. Stathopoulos, K. Wu; SIAM J. Sci. Comput. 23, 2165 (2002)
  //
  //
  if (_nConverged == _nEigen) {
    //if converged, S is empty.
    //Cannot orthonormalize empty matrix.
    return;
  }
  bool orthogonal = false;
  int nIter = 0;
  do {

    //calculate offDiagonalSum as measure for orthogonality
    //ToDo: There must be a more elegant way to measure orthogonality
    Eigen::MatrixXd thisShouldBeIdentity = S.transpose() * KS;
    double offDiagonalSum = 0;
    for (int i = 0; i < thisShouldBeIdentity.rows(); ++i) {
      for (int j = 0; j < thisShouldBeIdentity.cols(); ++j) {
        if (i != j) {
          offDiagonalSum += std::abs(thisShouldBeIdentity(i,j));
        }
      }
    }
    if (offDiagonalSum < 1.0e-12) {
      orthogonal = true;
//      print((std::string) "\n K-orthogonality ensured after " + nIter + " iterations" );
    }
    nIter += 1;
    //SVQB orthonormalization.
    Eigen::MatrixXd SVQB = S.transpose()  * KS;
    //
    VectorXd diagonalOfSVQB = SVQB.diagonal();
    for (int i = 0; i < diagonalOfSVQB.rows(); ++i) {
      diagonalOfSVQB(i) = 1 / std::sqrt(diagonalOfSVQB(i));
    }
    SVQB = diagonalOfSVQB.asDiagonal() * SVQB * diagonalOfSVQB.asDiagonal();
    SelfAdjointEigenSolver<MatrixXd> eigenvalueSolver(SVQB.cols());
    eigenvalueSolver.compute(SVQB);
    VectorXd eigenvaluesOfSVQB = eigenvalueSolver.eigenvalues();
    Eigen::MatrixXd eigenvectorsOfSVQB = eigenvalueSolver.eigenvectors();
    for (int i = 0; i < eigenvaluesOfSVQB.rows(); ++i) {
      if (eigenvaluesOfSVQB(i) < 1.0e-4) {
        eigenvaluesOfSVQB(i) = 1.0e-4;
      }
      eigenvaluesOfSVQB(i) = 1 / std::sqrt(eigenvaluesOfSVQB(i));
    }
    S = S * diagonalOfSVQB.asDiagonal() * eigenvectorsOfSVQB * eigenvaluesOfSVQB.asDiagonal();
    KS = KS * diagonalOfSVQB.asDiagonal() * eigenvectorsOfSVQB * eigenvaluesOfSVQB.asDiagonal();
    MKS = MKS * diagonalOfSVQB.asDiagonal() * eigenvectorsOfSVQB * eigenvaluesOfSVQB.asDiagonal();
    if ((nIter > 11 and offDiagonalSum > 1.0e-3) or
        (nIter >= 100)) {
      JacobiSVD<MatrixXd> svd(SVQB);
      double cond = svd.singularValues()(0)
          / svd.singularValues()(svd.singularValues().size()-1);
      print((std::string)" condition number = " + cond);
      throw SerenityError("could not k-orthonormalize W. It seems that W is ill conditioned.");
    }
  } while (orthogonal == false);
}

void KDavidsonEigenvalueSolver::subspaceCollapse() {
  //Will collapse the subspace dimension to _nEigen.
  //Not needed if there is enough memory available,
  //will slow convergence significantly.
  if (_subspaceDimension >= _maxSubspaceDimension) {
    print((std::string) "perform subspace collapse");
    _guessVectors[0] = _eigenvectors[0];
    _sigmaVectors[0] = _eigenvectors[1];
    _sigmaVectors[1] = _YM;
    _subspaceDimension = _nEigen;
    //Current expansion vectors are needed for toHDF5 of IterativeEigenvalueSolver
    _subspaceMatrix = _sigmaVectors[1].transpose() * _sigmaVectors[0];
    //Compute eigenpairs of reduced eigenvalue problem.
    //Use qr eigenvalue solver of eigen library.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(_subspaceDimension);
    eigenvalueSolver.compute(_subspaceMatrix);
    _subspaceEigenvalues = eigenvalueSolver.eigenvalues();
    _subspaceEigenvectors = eigenvalueSolver.eigenvectors();
  }
}

void KDavidsonEigenvalueSolver::iterate() {

  //
  //line 4: Form subspace matrix and solve reduced eigenvalue problem
  //
  _subspaceMatrix = _sigmaVectors[1].transpose() * _sigmaVectors[0];
  //Compute eigenpairs of reduced eigenvalue problem.
  //Use qr eigenvalue solver of eigen library.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(_subspaceDimension);
  eigenvalueSolver.compute(_subspaceMatrix);
  _subspaceEigenvalues = eigenvalueSolver.eigenvalues();
  _subspaceEigenvectors = eigenvalueSolver.eigenvectors();

  //
  //line 5: Get eigenvalues and compute trial vectors
  //
  _eigenvalues.resize(_nEigen);
  _eigenvalues.array() = _subspaceEigenvalues.block(0,0,_nEigen,1).array().sqrt();
  //compute trial vectors
  computeTrialVectors(_eigenvectors[0], _subspaceEigenvectors,_guessVectors[0]);
  computeTrialVectors(_eigenvectors[1],_subspaceEigenvectors,_sigmaVectors[0]);
  computeTrialVectors(_YM,_subspaceEigenvectors,_sigmaVectors[1]);
  //
  //lines 6 - 13:  calculate residual and correction vectors
  //
  computeResidualVectors(
      _residualVectors[0],
      _sigmaVectors[1],
      _guessVectors[0],
      _subspaceEigenvectors,
      _subspaceEigenvalues);
  checkForConvergence(_residualVectors);
  computeCorrectionVectors(_residualVectors[0],_W,_diagonalElements);

  //
  //lines 14 - 16: perform subspace collapse
  //
  subspaceCollapse();
  //
  //line 17: calculate new sigma vectors
  //
  _W = _W - _guessVectors[0] * ( _sigmaVectors[0].transpose() * _W);
  Eigen::MatrixXd wk = _KSigmaCalculator(_W);
  Eigen::MatrixXd wmk = _MSigmaCalculator(wk);

  //
  //line 18: k-Orthogonalize (iteratively to ensure k-orthogonality)
  //
  kOrthogonalize(_W,wk,wmk);

  //
  //line 19: Update guess vectors
  //
  int oldDimension = _subspaceDimension;
  int nAppend = _W.cols();
  _subspaceDimension = oldDimension + nAppend;
  _guessVectors[0].conservativeResize(_nDimension,_subspaceDimension);
  _guessVectors[0].block(0,oldDimension,_nDimension,nAppend) = _W;
  _sigmaVectors[0].conservativeResize(_nDimension,_subspaceDimension);
  _sigmaVectors[0].block(0,oldDimension,_nDimension,nAppend) = wk;
  _sigmaVectors[1].conservativeResize(_nDimension,_subspaceDimension);
  _sigmaVectors[1].block(0,oldDimension,_nDimension,nAppend) = wmk;

  //print some data
  printIterationData(_residualVectors,_eigenvalues);
}

void KDavidsonEigenvalueSolver::postProcessing() {
  //line 21:
  //calculate other component of the eigenvector
  VectorXd inverseEigenvalues(_nEigen);
  for (unsigned int i = 0; i < _nEigen; ++i) {
    if (_eigenvalues(i) < 1.0e-6 ) {
      std::cout << "Warning: Found negative eigenvalue?!";
    }
    inverseEigenvalues(i) = 1.0 / _eigenvalues(i);
  }
  _eigenvectors[1] = _eigenvectors[1] * inverseEigenvalues.asDiagonal();
  //Remove nan's for analysis
  auto nans = _eigenvalues.array().isNaN();
  for (unsigned int i = 0; i < _nEigen; ++i) {
    if (nans(i) == 1) {
      _eigenvectors[0].col(i).setZero();
      _eigenvectors[1].col(i).setZero();
    }
  }

  //get u and v (according to article by Vecharynski et. al), i.e. obtain
  //X and Y component of the eigenvector  (Y,X)^T
  Eigen::MatrixXd oldX = _eigenvectors[0];
  Eigen::MatrixXd oldY = _eigenvectors[1];
  double factor = 1.0 / std::sqrt(2);
  _eigenvectors[0] =  factor*(oldY + oldX);
  _eigenvectors[1] =  factor*(oldY - oldX);

  auto& x = _eigenvectors[0];
  auto& y = _eigenvectors[1];
  Eigen::MatrixXd tmp = x.transpose()*x - y.transpose()*y;
  Eigen::VectorXd norms = tmp.colwise().norm();
  for (unsigned int i = 0; i < norms.rows(); ++i) {
    norms(i) = 1.0 / std::sqrt(norms(i));
  }
  x *= norms.asDiagonal();
  y *= norms.asDiagonal();

  _hasBeenRun = true;
}


Eigen::MatrixXd& KDavidsonEigenvalueSolver::getX() {
  if (_hasBeenRun == false) {
    solve();
  }
  return _eigenvectors[0];
}

Eigen::MatrixXd& KDavidsonEigenvalueSolver::getY() {
  if (_hasBeenRun == false) {
    solve();
  }
  return _eigenvectors[1];
}

void KDavidsonEigenvalueSolver::initializeKDavidson() {
  //print some data
  printHeader("K-Davidson eigenvalue solver");
  print((std::string)"Maximum subspace dimension          : " + _maxSubspaceDimension);
  //What is optimized here are some eigenvectors of (A+B)(A-B).
  //Since the diagonal elements are orbital energy differences here,
  //and if we assume that B_ii is much smaller than A_ii, we can in a good
  //approximation express the diagonal elements by squared orbital energy
  //differences.
  _diagonalElements = _diagonalElements.array().square();
  //Set dimensions
  _guessVectors.resize(1);
  _eigenvectors.resize(2);
  _residualVectors.resize(1);
  _sigmaVectors.resize(2);
  _method.resize(1);
  _eigenvectors[0].resize(_nDimension, _nEigen);
  _eigenvectors[1].resize(_nDimension, _nEigen);
  _YM.resize(_nDimension, _nEigen);

  //Try to initialize from file, else normal start
  if (!initializeFromHDF5()) {
    _guessVectors[0].resize(_nDimension,_nEigen);
    // line 1 of Algorithm 5: Produce k-orthogonal set of guess vectors
    computeInitialGuessVectors(_guessVectors[0], _diagonalElements);
    //This expensive steps are necessary to k-Orthogonalize the initial guess
    _sigmaVectors[0] = _KSigmaCalculator(_guessVectors[0]);
    _sigmaVectors[1] = _MSigmaCalculator(_sigmaVectors[0]);
    //kOrthogonalize.
    kOrthogonalize(_guessVectors[0], _sigmaVectors[0], _sigmaVectors[1]);
  } else {
    //If more eigenvalues than stored are requested, we need to add additional guess
    //and sigma vectors
    if (_nEigen > _guessVectors[0].cols()) {
      //ToDo: Reuse guess vectors from HDF5 file...
      _guessVectors[0].resize(_nDimension,_nEigen);
      // line 1 of Algorithm 5: Produce k-orthogonal set of guess vectors
      computeInitialGuessVectors(_guessVectors[0], _diagonalElements);
      //This expensive steps are necessary to k-Orthogonalize the initial guess
      _sigmaVectors[0] = _KSigmaCalculator(_guessVectors[0]);
      _sigmaVectors[1] = _MSigmaCalculator(_sigmaVectors[0]);
      //kOrthogonalize.
      kOrthogonalize(_guessVectors[0], _sigmaVectors[0], _sigmaVectors[1]);
    }
  }
  _subspaceDimension = _guessVectors[0].cols();
}


} /* namespace Serenity */
