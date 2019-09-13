/**
 * @file EigenvalueSolver.cpp
 *
 * @date May 30, 2019
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
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"

/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "parameters/Constants.h"
#include "misc/SerenityError.h"
#include <fstream>
namespace Serenity {

  EigenvalueSolver::EigenvalueSolver(
    bool printResponseMatrix,
    unsigned nDimension,
    unsigned nEigen,
    Eigen::VectorXd& diagonal,
    double convergenceCriterion,
    unsigned maxIterations,
    unsigned maxSubspaceDimension,
    unsigned initialSubspace,
    Options::RESPONSE_PROBLEM responseType,
    std::function<std::unique_ptr<std::vector<Eigen::MatrixXd> >(
      std::vector<Eigen::MatrixXd>& guessVectors)> sigmaCalculator,
    std::shared_ptr<std::vector<Eigen::MatrixXd> > initialGuess) : 
    IterativeSolver(
      responseType == Options::RESPONSE_PROBLEM::RPA ? 2 : 1,
      nDimension,nEigen,diagonal,convergenceCriterion,
      maxIterations,maxSubspaceDimension,sigmaCalculator,initialGuess),
    _printResponseMatrix(printResponseMatrix),
    _initialSubspace(initialSubspace),
    _responseType(responseType){
      this->initialize();
  }

  void EigenvalueSolver::initialize(){

    printHeader("Eigenvalue Solver");
    if(_responseType == Options::RESPONSE_PROBLEM::TDA){
      printf("\n    Solving Hermitian CIS/TDA problem.\n");
    }else if(_responseType == Options::RESPONSE_PROBLEM::TDDFT){
      printf("\n    Solving Hermitian TDDFT problem.\n");
    }else{
      printf("\n    Solving non-Hermitian RPA problem.\n");
    }

    //Set dimensions
    _sigmaVectors.resize(_nSets);
    _guessVectors.resize(_nSets);
    _eigenvectors.resize(_nSets);
    _residualVectors.resize(_nSets);
    _expansionVectors.resize(_nSets);
    _correctionVectors.resize(_nSets);
    _residualNorms = Eigen::VectorXd::Zero(_nEigen);
    _eigenvalues = Eigen::VectorXd::Zero(_nEigen);

    //Initialize all matrices
    for(unsigned iSet = 0; iSet < _nSets; ++iSet){
      _guessVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension,_initialSubspace);
      _eigenvectors[iSet] = Eigen::MatrixXd::Zero(_nDimension,_nEigen);
      _residualVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension,_nEigen);
      _correctionVectors[iSet] = Eigen::MatrixXd::Zero(_nDimension,_nEigen);
    }

    //Set guess vectors to 
    // - the initial guess if one was passed to the eigenvalue solver
    // - unit vectors in case of fresh calculation
    if(_initialGuess){
      if(_nSets > (*_initialGuess).size()) (*_initialGuess).resize(_nSets);
      for(unsigned iSet = 0; iSet < _nSets; ++iSet){
        _guessVectors[iSet] = (*_initialGuess)[iSet];
      }
      Orthogonalization::modifiedGramSchmidtLinDep(_guessVectors,0.0);
    }else{
      std::vector<std::pair<double,unsigned> > diagPair;
      for(unsigned ia = 0; ia < _nDimension; ++ia){
        diagPair.push_back(std::pair<double,unsigned>(_diagonal(ia), ia));
      }
      std::stable_sort(diagPair.begin(), diagPair.end());
      for(unsigned iSet = 0; iSet < _nSets; ++iSet){
        for(unsigned i = 0; i < _initialSubspace; ++i){
          _guessVectors[iSet](diagPair[i].second, i) = 1.0;
        }
      }
    }

    //Get initial sigma vectors
    _sigmaVectors = (*_sigmaCalculator(_guessVectors));
  } /* this->initialize() */

  void EigenvalueSolver::iterate(){
    
    //Set current subspace size
    _subDim = _guessVectors[0].cols();
    
    //Solve subspace problem
    //At the end of this section, the matrix _subspaceEigenvectors will contain the current expansion
    //coefficients of this iteration sorted according to their respective eigenvalue
    if(_nSets == 1){
      //Construct Hermitian subspace problem and solve it exactly
      Eigen::MatrixXd odiag = (_guessVectors[0].transpose() * _guessVectors[0]).diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
      _subspaceMatrix = odiag * _guessVectors[0].transpose() * _sigmaVectors[0] * odiag;
      if(_printResponseMatrix) {
        std::ofstream file("ResponseMatrix_"+std::to_string(this->_nIter)+".txt");
        if (file.is_open())
        {
          file << _subspaceMatrix;
        }
        file.close();
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(_subspaceMatrix);
      _subspaceEigenvectors = odiag * subspaceSolver.eigenvectors();
      _subspaceEigenvalues = subspaceSolver.eigenvalues();
    }else{
      //Construct symplectic subspace problem and solve it exactly
      _subspaceMatrix = Eigen::MatrixXd::Zero(2*_subDim,2*_subDim);
      _subspaceMatrix.topLeftCorner(_subDim,_subDim) = _guessVectors[0].transpose() * _sigmaVectors[0];
      _subspaceMatrix.bottomRightCorner(_subDim,_subDim) = _guessVectors[1].transpose() * _sigmaVectors[1];

      Eigen::MatrixXd overlap = Eigen::MatrixXd::Zero(2*_subDim,2*_subDim);
      overlap.topRightCorner(_subDim,_subDim) = _guessVectors[0].transpose() * _guessVectors[1];
      overlap.bottomLeftCorner(_subDim,_subDim) = _guessVectors[1].transpose() * _guessVectors[0];

      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> subspaceSolver(_subspaceMatrix,overlap,true);
      Eigen::MatrixXd unsortedEigenvectors = subspaceSolver.eigenvectors().real();
      Eigen::VectorXd unsortedEigenvalues = subspaceSolver.eigenvalues().real();

      //Sort (only positive) eigenvalues in ascending order
      std::vector<std::pair<double,unsigned> > evPair(0);
      for(unsigned i = 0; i < 2 * _subDim; ++i){
        if (unsortedEigenvalues(i) > 0) evPair.push_back(std::make_pair(unsortedEigenvalues(i),i));
      }
      std::stable_sort(evPair.begin(), evPair.end());

      _subspaceEigenvectors.resize(2*_subDim,_subDim);
      _subspaceEigenvalues.resize(_subDim);
      //Reassemble eigenvectors
      for(unsigned i = 0; i < _subDim; ++i){
        _subspaceEigenvalues(i) = evPair[i].first;
        _subspaceEigenvectors.col(i) = unsortedEigenvectors.col(evPair[i].second);
      }
    }

    //If there was an initial guess, the eigenvalue solver will try to follow _nEigen roots
    //that were given to it based on the overlap of the ritz vectors with the initial guess
    if(_initialGuess){
      Eigen::MatrixXd& rootsToFollow = _nIter == 1 ? (*_initialGuess)[0] : _eigenvectors[0];
      Eigen::MatrixXd overlap = (_guessVectors[0] * _subspaceEigenvectors.topRows(_subDim)).transpose() * rootsToFollow;
      for(unsigned i = 0; i < _nEigen; ++i){
        unsigned iMax;
        overlap.col(i).cwiseAbs().maxCoeff(&iMax);
        //Swap eigenvalues and expansion vectors
        _subspaceEigenvalues.row(i).swap(_subspaceEigenvalues.row(iMax));
        _subspaceEigenvectors.col(i).swap(_subspaceEigenvectors.col(iMax));
      }
    }
 
    //Sort eigenvalues in ascending order that were found to have the best overlap
    std::vector<std::pair<double,unsigned> > evPair(0);
    for(unsigned i = 0; i < _nEigen; ++i){
      evPair.push_back(std::make_pair(_subspaceEigenvalues(i),i));
    }
    std::stable_sort(evPair.begin(), evPair.end());

    //Reassemble eigenvalues and expansion coefficients
    for(unsigned iSet = 0; iSet < _nSets; ++iSet){
      _expansionVectors[iSet].resize(_subDim,_nEigen);
    }
    for(unsigned i = 0; i < _nEigen; ++i){
      _eigenvalues(i) = evPair[i].first;
      for(unsigned iSet = 0; iSet < _nSets; ++iSet){
        _expansionVectors[iSet].col(i) = _subspaceEigenvectors.middleRows(iSet*_subDim,_subDim).col(evPair[i].second);
      }
    }

    //Ritz vectors
    for(unsigned iSet = 0; iSet < _nSets; ++iSet){
      _eigenvectors[iSet] = _guessVectors[iSet] * _expansionVectors[iSet];
    }

    //Residual vectors and norms
    _residualNorms.setZero();
    for(unsigned iSet = 0; iSet < _nSets; ++iSet){
      _residualVectors[iSet] = _sigmaVectors[iSet] * _expansionVectors[iSet] 
        - _eigenvectors[_nSets == 1 ? 0 : iSet == 0 ? 1 : 0] * _eigenvalues.asDiagonal();
      _residualNorms += _residualVectors[iSet].colwise().norm() / std::sqrt(_responseType == Options::RESPONSE_PROBLEM::TDDFT ? 4 : _nSets);
    }

    //Check convergence
    _nConverged = 0;
    for(unsigned i = 0; i < _nEigen; ++i){
      if(_residualNorms(i) < _convergenceCriterion) ++_nConverged;
    }

    //Print iteration data
    if (_nIter == 1) printTableHead("\n   It.    Dimension    Converged    Max. Norm");
    printf("%5i %11i %12i %14.2e\n",_nIter,(int)_subDim,_nConverged,_residualNorms.maxCoeff());

    //Skip the rest of this iteration if converged or no further iterations are planned
    if(_nConverged == _nEigen || _nIter == _maxIterations) return;

    //Reset correction vectors
    for(unsigned iSet = 0; iSet < _nSets; ++iSet){
      _correctionVectors[iSet].resize(_nDimension,_nEigen-_nConverged);
      _correctionVectors[iSet].setZero();
    }
  
    //Use diagonal approximation to precondition residual vectors
    int iCount = -1;
    for(unsigned j = 0; j < _nEigen; ++j){
      if(_residualNorms(j) < _convergenceCriterion) continue;
      ++iCount;
      for(unsigned ia = 0; ia < _nDimension; ++ia){
        if(_responseType == Options::RESPONSE_PROBLEM::TDA){
          double aj = 1 / (_eigenvalues(j) - _diagonal(ia));
          _correctionVectors[0](ia,iCount) = aj * _residualVectors[0](ia,j);
        }else if(_responseType == Options::RESPONSE_PROBLEM::TDDFT){
          double aj = 1 / (_eigenvalues(j) - _diagonal(ia)*_diagonal(ia));
          _correctionVectors[0](ia,iCount) = aj * _residualVectors[0](ia,j);
        }else{
          double aj = 1 / (_diagonal(ia) - _eigenvalues(j) / _diagonal(ia) * _eigenvalues(j));
          double bj = aj  * _eigenvalues(j) / _diagonal(ia);
          _correctionVectors[0](ia,iCount) = aj * _residualVectors[0](ia,j) + bj * _residualVectors[1](ia,j);
          _correctionVectors[1](ia,iCount) = bj * _residualVectors[0](ia,j) + aj * _residualVectors[1](ia,j);
        }
      }
    }

    //Expand subspace
    this->expandSubspace();
  } /* this->iterate() */

  void EigenvalueSolver::postProcessing(){

    if(_responseType == Options::RESPONSE_PROBLEM::TDA){
      //Do nothing
    }else{
      //In case of TDDFT, obtain left eigenvectors
      if(_responseType == Options::RESPONSE_PROBLEM::TDDFT){
        printf(" Transforming for TDDFT.\n\n");
        //Obtain eigenvalues
        _eigenvalues = _eigenvalues.cwiseSqrt();
        //Obtain (X+Y)
        _eigenvectors[0] = _diagonal.cwiseSqrt().asDiagonal() * _eigenvectors[0];
        //Obtain (X-Y)
        _eigenvectors.push_back(_diagonal.cwiseInverse().asDiagonal() * _eigenvectors[0] * _eigenvalues.asDiagonal());
      }
      //Obtain X and Y from (X+Y) and (X-Y)
      Eigen::MatrixXd xpy = _eigenvectors[0];
      Eigen::MatrixXd xmy = _eigenvectors[1];
      _eigenvectors[0] = (xpy + xmy) / std::sqrt(2);
      _eigenvectors[1] = (xpy - xmy) / std::sqrt(2);

      //Normalize
      Eigen::MatrixXd overlap = _eigenvectors[0].transpose()*_eigenvectors[0] - _eigenvectors[1].transpose()*_eigenvectors[1];
      Eigen::VectorXd norms = overlap.colwise().norm();
      for(unsigned i = 0; i < norms.size(); ++i){
        norms(i) = 1 / std::sqrt(norms(i));
      }
      _eigenvectors[0] *= norms.asDiagonal();
      _eigenvectors[1] *= norms.asDiagonal();
    }

    //Print eigenvalues with norms
    printf("--------------------------------------------\n");
    printf("   root        eigenvalue       res. norm   \n");
    printf("           (a.u.)      (eV)                 \n");
    printf("--------------------------------------------\n");
    for(unsigned i = 0; i < _nEigen; ++i){
      printf("%6i %11.6f %9.4f %12.3e\n",
      i + 1,
      _eigenvalues(i),
      _eigenvalues(i) * HARTREE_TO_EV,
      _residualNorms(i));
    }
    printf("\n");
  } /* this->postProcessing() */
    
} /* namespace Serenity */
