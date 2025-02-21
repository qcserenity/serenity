/**
 * @file SQNM.cpp
 *
 * @date Oct 09, 2024
 * @author Kasibek Zumataev, reworked by Thorben Wiegmann
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
/* Include Serenity Internal Headers */
#include "math/optimizer/SQNM.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

SQNM::SQNM(Eigen::VectorXd& parameters, unsigned int historyLength, double epsilon, double alpha,
           double energyThreshold, double trustRadius)
  : Optimizer(parameters),
    _historyLength(historyLength),
    _epsilon(epsilon),
    _alpha(alpha),
    _alphaStart(alpha),
    _eThresh(energyThreshold),
    _trustRadius(trustRadius) {
}

Eigen::VectorXd SQNM::scaleStep(Eigen::VectorXd step) {
  if (step.maxCoeff() > _trustRadius) {
    OutputControl::mOut << "Trust radius exceeded. Scaling step accordingly..." << std::endl;
    step *= _trustRadius / step.maxCoeff();
  }
  return step;
}

std::vector<Eigen::VectorXd> SQNM::getDifferences(std::vector<Eigen::VectorXd> vec) {
  std::vector<Eigen::VectorXd> resVec;
  for (unsigned int i = 0; i < vec.size() - 1; i++) {
    resVec.push_back(vec[i] - vec[i + 1]);
  }
  if (resVec.size() != vec.size() - 1) {
    throw SerenityError("Calculation of difference vector returned a wrong dimensionality!");
  }
  return resVec;
}

Eigen::MatrixXd SQNM::calcNorms(std::vector<Eigen::VectorXd> vec) {
  Eigen::MatrixXd resMat = Eigen::MatrixXd(this->_parameters.size(), vec.size()).setZero();
  for (unsigned int i = 0; i < vec.size(); i++) {
    resMat.col(i) = vec[i] / vec[i].norm();
  }
  return resMat;
}

void SQNM::optimize(std::function<bool(const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                       std::shared_ptr<Eigen::MatrixXd> hessian, bool print)>
                        updateFunction,
                    std::shared_ptr<unsigned int> nRejected) {
  // initialize variables
  double value = std::numeric_limits<double>::infinity();
  double oldValue = std::numeric_limits<double>::infinity();
  Eigen::VectorXd gradients = Eigen::VectorXd(this->_parameters.size()).setZero();
  Eigen::VectorXd oldGradients = Eigen::VectorXd(this->_parameters.size()).setZero();
  Eigen::VectorXd oldPrecon = Eigen::VectorXd(this->_parameters.size()).setZero();
  bool converged = false;
  // check convergence once before the loop starts and calculate first energy and gradients
  converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
  // check system and gradients size
  if (this->_parameters.size() != gradients.size()) {
    throw SerenityError("Gradients do not match the system!");
  }
  // start while loop of optimization
  while (true) {
    // check convergence in the beginning
    if (converged) {
      break;
    }
    // append current coordinates and gradients
    _coordList.push_back(this->_parameters);
    _gradientList.push_back(gradients);
    // perform steepest descent step, if history is too short for SQNM step
    if (_coordList.size() < 3) {
      Eigen::VectorXd step = _alpha * gradients;
      this->_parameters -= this->scaleStep(step);
      // check convergence
      converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
    }
    // perform SQNM step
    else {
      // get displacements, normalized displacements, and gradient differences
      auto displacements = this->getDifferences(_coordList);
      auto normedDisplacements = this->calcNorms(displacements);
      auto gradDiffs = this->getDifferences(_gradientList);
      unsigned int nHist = normedDisplacements.cols();
      // calculate overlap matrix, eigenvalues and eigenvectors
      Eigen::MatrixXd ovlp = normedDisplacements.transpose() * normedDisplacements;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ovlpSolver(ovlp);
      Eigen::VectorXd ovlpValues = ovlpSolver.eigenvalues();
      Eigen::MatrixXd ovlpVectors = ovlpSolver.eigenvectors();
      // determine significant eigenvalues
      unsigned int cutoff = 0;
      for (unsigned int i = 0; i < ovlpValues.rows(); i++) {
        // eigenvalues are sorted from small to large by Eigen, thus the first significant eigenvalue will be followed
        // by other significant eigenvalues
        if (ovlpValues(i) / ovlpValues.maxCoeff() > _epsilon) {
          cutoff = i;
          break;
        }
      }
      // truncate eigenvalues and eigenvectors to only contain significant contributions
      Eigen::VectorXd sigValues = ovlpValues.tail(ovlp.rows() - cutoff);
      Eigen::MatrixXd sigVectors = ovlpVectors.rightCols(ovlpVectors.cols() - cutoff);
      unsigned int nDim = sigValues.rows();
      // define significant subspace
      Eigen::MatrixXd sigSubspace = normedDisplacements * sigVectors;
      for (unsigned int i = 0; i < nDim; i++) {
        sigSubspace.col(i) *= 1 / std::sqrt(sigValues(i));
      }
      // calculate gradient differences in the significant subspace
      Eigen::MatrixXd sigGradDiff = Eigen::MatrixXd(this->_parameters.size(), nDim).setZero();
      for (unsigned int i = 0; i < nDim; i++) {
        for (unsigned int k = 0; k < nHist; k++) {
          sigGradDiff.col(i) += (sigVectors(k, i) / displacements[k].norm()) * gradDiffs[k];
        }
        sigGradDiff.col(i) *= 1 / std::sqrt(sigValues(i));
      }
      // calculate approximate Hessian and symmetrize it, calculate eigenvalues and eigenvectors
      Eigen::MatrixXd hessian = sigGradDiff.transpose() * sigSubspace;
      hessian = 0.5 * (hessian + hessian.transpose()).eval();
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessianSolver(hessian);
      Eigen::VectorXd hessianValues = hessianSolver.eigenvalues();
      Eigen::MatrixXd hessianVectors = hessianSolver.eigenvectors();
      // calculate search directions in the subspace
      Eigen::MatrixXd searchDirections = sigSubspace * hessianVectors;
      // calculate projection matrix and orthogonal gradient part
      Eigen::MatrixXd projector = Eigen::MatrixXd(this->_parameters.size(), this->_parameters.size()).setZero();
      for (unsigned int i = 0; i < nDim; i++) {
        projector += searchDirections.col(i) * searchDirections.col(i).transpose();
      }
      Eigen::MatrixXd identity = Eigen::MatrixXd(this->_parameters.size(), this->_parameters.size()).setIdentity();
      Eigen::VectorXd orthoGrad = (identity - projector) * gradients;
      // calculate residues and modified eigenvalues of the hessian
      Eigen::VectorXd residues = Eigen::VectorXd(nDim).setZero();
      Eigen::MatrixXd preMat = sigGradDiff * hessianVectors;
      for (unsigned int i = 0; i < nDim; i++) {
        residues(i) = (preMat.col(i) - hessianValues(i) * searchDirections.col(i)).norm();
      }
      Eigen::VectorXd modifiedValues = Eigen::VectorXd(nDim).setZero();
      for (unsigned int i = 0; i < nDim; i++) {
        modifiedValues(i) = std::sqrt(std::pow(hessianValues(i), 2) + std::pow(residues(i), 2));
      }
      // calculate preconditioned gradient of the subspace
      Eigen::VectorXd preconGrad = Eigen::VectorXd(this->_parameters.size()).setZero();
      for (unsigned int i = 0; i < nDim; i++) {
        preconGrad += (gradients.dot(searchDirections.col(i)) / modifiedValues(i)) * searchDirections.col(i);
      }
      // adjust alpha
      double scaling = gradients.dot(preconGrad) / (gradients.norm() * preconGrad.norm());
      if (scaling > 0.2) {
        _alpha *= 1.1;
      }
      else {
        _alpha *= 0.85;
      }
      // calculate total preconditioned gradient
      Eigen::VectorXd totalGrad = preconGrad + _alpha * orthoGrad;
      this->_parameters -= this->scaleStep(totalGrad);
      oldValue = value;
      oldGradients = gradients;
      converged = updateFunction(this->_parameters, value, gradients, nullptr, true);
      // check for fallback
      if (value > (oldValue + _eThresh) && _alpha > (_alphaStart / 10.0)) {
        // step is not accepted
        _coordList.clear();
        _gradientList.clear();
        _alpha /= 2.0;
        value = oldValue;
        gradients = oldGradients;
        OutputControl::mOut << "Step not accepted. Resetting coordinates to previous step." << std::endl;
        *(nRejected) = *(nRejected) + 1;
        // reverse step taken by SQNM - check for trust radius
        this->_parameters += this->scaleStep(totalGrad);
        converged = false;
      }
    } // SQNM algorithm
    if (_coordList.size() > _historyLength) {
      while (_coordList.size() > _historyLength) {
        _coordList.erase(_coordList.begin());
      }
    }
    if (_gradientList.size() > _historyLength) {
      while (_gradientList.size() > _historyLength) {
        _gradientList.erase(_gradientList.begin());
      }
    }
  } // while loop
} // optimize function
} /* namespace Serenity */
