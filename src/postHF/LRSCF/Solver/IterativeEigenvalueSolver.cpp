/**
 * @file IterativeEigenvalueSolver.cpp
 *
 * @date Oct 17, 2016
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
#include "postHF/LRSCF/Solver/IterativeEigenvalueSolver.h"
/* Include Serenity Internal Headers */
#include "io/HDF5.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <iostream>


namespace Serenity {

  IterativeEigenvalueSolver::IterativeEigenvalueSolver(
      unsigned int nDimension,
      unsigned int nEigen,
      Eigen::VectorXd& diagonalElements,
      double convergenceCriterion,
      unsigned int nMaxIterations,
      std::string loadPath,
      std::string filePath,
      std::string id):
        _nDimension(nDimension),
        _nEigen(nEigen),
        _diagonalElements(diagonalElements),
        _convergenceCriterion(convergenceCriterion),
        _nMaxIterations(nMaxIterations),
        _loadPath(loadPath),
        _filePath(filePath),
        _id(id),
        _restart(false){
    //Set dimension of objects to zero. The dimension has to be set in actual implementations
    //of IterativeEigenvalueSolver
    _guessVectors.resize(0);
    _sigmaVectors.resize(0);
    _residualVectors.resize(0);
    _eigenvectors.resize(0);
    if (_filePath != "") _restart = true;
  }

  std::vector<Eigen::MatrixXd >& IterativeEigenvalueSolver::getEigenvectors() {
    //calculate eigenvectors if not already done
    if (_hasBeenRun == false) {
      solve();
    }
    return _eigenvectors;
  }

  Eigen::VectorXd& IterativeEigenvalueSolver::getEigenvalues() {
    //calculate eigenvalues if not already done
    if (_hasBeenRun == false) {
      solve();
    }
    return _eigenvalues;
  }

  std::vector<Eigen::MatrixXd >& IterativeEigenvalueSolver::getResiduals() {
    //calculate if not already done
    if (_hasBeenRun == false) {
      solve();
    }
    return _residualVectors;
  }

  void IterativeEigenvalueSolver::solve() {
//    Timings::takeTime("Iterative eigenvalue solver");
    do {
      _nIter += 1;
      iterate();
      //Write guess and sigma vectors to temporary file for restart
      if (_restart) toHDF5();
      //exit if maximum number of iterations is reached
      if ( _nIter >= _nMaxIterations) {
        print((std::string)"\n Reached maximum number of iterations \n");
        break;
      }
      if (_nConverged == _nEigen) {
        _converged = true;
      }
    } while (_converged == false);
    if (_converged) {
      print((std::string)"\n Reached convergence criterion \n");
    } else {
      WarningTracker::printWarning("Warning: convergence criterion not reached",true);
    }
    postProcessing();
    _hasBeenRun = true;
//    Timings::timeTaken("Iterative eigenvalue solver");
  }

  bool& IterativeEigenvalueSolver::isConverged() {
    if (_hasBeenRun == false) {
      solve();
    }
    return _converged;
  }

  void IterativeEigenvalueSolver::computeInitialGuessVectors(
      Eigen::MatrixXd& guessVectors,
      Eigen::VectorXd& diagonalElements) {
    guessVectors = Eigen::MatrixXd::Zero(_nDimension, _nEigen);
    //write diagonal elements to pair vector to get lowest elements and corresponding indices
    std::vector<std::pair<double, int>> pairVector;
    for (unsigned int i = 0; i < _nDimension; ++i) {
      pairVector.push_back(std::pair<double, int>(diagonalElements(i), i));
    }
    std::stable_sort(pairVector.begin(), pairVector.end());
    for (unsigned int i = 0; i < _nEigen; ++i) {
      guessVectors(pairVector[i].second, i) = 1.0;
    }
  }

  void IterativeEigenvalueSolver::computeTrialVectors(
    Eigen::MatrixXd& trialVectors,
    Eigen::MatrixXd& expansionCoefficients,
    Eigen::MatrixXd& guessVectors) {
    //This is what actually is done here is
    //
    //  trialVectors = guessVectors * expansionCoefficients;
    //
    //However, this would be slower in most cases since we only need _nEigen columns.
    trialVectors = Eigen::MatrixXd::Zero(_nDimension,_nEigen);
    for (unsigned int i = 0; i < _nEigen; ++i) {
      for (unsigned int j = 0; j < expansionCoefficients.cols(); ++j) {
        trialVectors.col(i) += expansionCoefficients(j, i) * guessVectors.col(j);
      }
    }
  }

  void IterativeEigenvalueSolver::computeResidualVectors(
      Eigen::MatrixXd& residualVectors,
      Eigen::MatrixXd& sigmaVectors,
      Eigen::MatrixXd& guessVectors,
      Eigen::MatrixXd& expansionCoefficients,
      Eigen::VectorXd& subspaceEigenvalues) {
    residualVectors = Eigen::MatrixXd::Zero(_nDimension,_nEigen);
    for (unsigned int i = 0; i < _nEigen; ++i) {
      for (unsigned int j = 0; j < subspaceEigenvalues.rows(); ++j) {
        residualVectors.col(i) += expansionCoefficients(j, i)
            * (sigmaVectors.col(j) - subspaceEigenvalues(i) * guessVectors.col(j));
      }
    }
  }

  void IterativeEigenvalueSolver::checkForConvergence(std::vector<Eigen::MatrixXd >& residualVectors) {
    _nConverged = 0;
    for (unsigned int i = 0; i < _nEigen; ++i) {
      double maxNorm = 0.0;
      for (unsigned int j = 0; j < residualVectors.size(); ++j) {
        maxNorm = std::max(maxNorm,residualVectors[j].col(i).norm());
      }
      if (maxNorm < _convergenceCriterion) {
        _nConverged += 1;
      }
    }
  }

  void IterativeEigenvalueSolver::printIterationData(
      std::vector<Eigen::MatrixXd >& residualVectors,
      Eigen::VectorXd& eigenvalues){
    std::cout << std::fixed;
    std::cout << std::setprecision(12);
    //Print table head in first iteration
    if (_nIter == 1) {
        print(" ");
        printTableHead(" It.    Dimension  #converged    max. residual norm");
    }
    (void)eigenvalues; //no warning, please , TODO print this argument?
      //calculate maxNorm of all residual vectors
      double maxNorm = 0.0;
      for (unsigned int i = 0; i < residualVectors.size(); ++i) {
        maxNorm = std::max(maxNorm,residualVectors[i].colwise().norm().maxCoeff());
      }
      std::cout << "    " << std::setw(3) << _nIter << "       "
          << std::setw(3) << _subspaceEigenvalues.rows() << "       " << std::setw(3)
          << _nConverged << "        " << std::setw(15) << maxNorm << std::endl;
  }

  void IterativeEigenvalueSolver::printHeader(std::string caption) {
      //print caption and some general information
      printBigCaption(caption);
      print((std::string)"Dimension of the eigenvalue problem : " + _nDimension);
      print((std::string)"Convergence threshold               : " + _convergenceCriterion);
      print((std::string)"Maximum number of iterations        : " + _nMaxIterations);
  }

  void IterativeEigenvalueSolver::toHDF5() {
    std::string fName = _filePath+".solver.h5";
    HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);
    for (unsigned int iGuess = 0; iGuess < _guessVectors.size(); ++iGuess) {
      std::string dataset = "guess"+std::to_string(iGuess);
      //Expand guessVectors, so that only _nEigen columns need to be stored.
      //This corresponds to a subspace collapse for the restarted calculation.
      Eigen::MatrixXd eigenvectors(_nDimension,_nEigen);
      if (_guessVectors[iGuess].cols() == _nEigen) {
        eigenvectors = _guessVectors[iGuess];
      } else  {
        computeTrialVectors(eigenvectors,_subspaceEigenvectors,_guessVectors[iGuess]);
      }
      HDF5::save(file,dataset.c_str(),eigenvectors);
    }
    for (unsigned int iSigma = 0; iSigma < _sigmaVectors.size(); ++iSigma) {
      std::string dataset = "sigma"+std::to_string(iSigma);
      //Expand sigmaVectors, so that only _nEigen columns need to be stored
      //This corresponds to a subspace collapse for the restarted calculation.
      Eigen::MatrixXd expandedSigmaVectors(_nDimension,_nEigen);
      if (_sigmaVectors[iSigma].cols() == _nEigen) {
        expandedSigmaVectors = _sigmaVectors[iSigma];
      } else {
        computeTrialVectors(expandedSigmaVectors,_subspaceEigenvectors,_sigmaVectors[iSigma]);
      }
      HDF5::save(file,dataset.c_str(),expandedSigmaVectors);
    }
    HDF5::save_scalar_attribute(file,"ID",_id);
    file.close();
  }

  bool IterativeEigenvalueSolver::initializeFromHDF5() {
    if (_sigmaVectors.size() == 0 || _guessVectors.size() == 0) {
      throw SerenityError("IterativeEigenvalueSolver has not been initialized correctly.");
    }
    try {
      std::string fName = _loadPath+".solver.h5";
      HDF5::Filepath name(fName);
      HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      HDF5::attribute_exists(file,"ID");
      HDF5::check_attribute(file,"ID",_id);
      for (unsigned int iGuess = 0; iGuess < _guessVectors.size(); ++iGuess) {
        std::string dataset = "guess"+std::to_string(iGuess);
        HDF5::dataset_exists(file,dataset.c_str());
        HDF5::load(file,dataset.c_str(),_guessVectors[iGuess]);
      }
      for (unsigned int iSigma = 0; iSigma < _sigmaVectors.size(); ++iSigma) {
        std::string dataset = "sigma"+std::to_string(iSigma);
        HDF5::dataset_exists(file,dataset.c_str());
        HDF5::load(file,dataset.c_str(),_sigmaVectors[iSigma]);
      }
      file.close();
      return true;
    } catch (...) {
      return false;
    }
  }
}
