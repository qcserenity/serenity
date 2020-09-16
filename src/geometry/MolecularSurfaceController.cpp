/**
 * @file MolecularSurfaceController.cpp
 *
 * @author Moritz Bensberg
 * @date May 25, 2020
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
#include "geometry/MolecularSurfaceController.h"
#include "grid/GridController.h"
#include "io/FormattedOutputStream.h"
#include "math/linearAlgebra/MatrixFunctions.h" //Pseudo inversion.

namespace Serenity {

MolecularSurfaceController::MolecularSurfaceController(std::shared_ptr<GridController> gridController,
                                                       std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string surfaceLabel)
  : _gridController(gridController), _normalVectors(std::move(normalVectors)), _surfaceLabel(surfaceLabel) {
  printInfo();
}

MolecularSurfaceController::MolecularSurfaceController(std::unique_ptr<Eigen::Matrix3Xd> coordinates,
                                                       std::unique_ptr<Eigen::VectorXd> weights,
                                                       std::unique_ptr<Eigen::Matrix3Xd> normalVectors, std::string surfaceLabel)
  : _gridController(
        std::make_shared<GridController>(std::unique_ptr<Grid>(new Grid(std::move(coordinates), std::move(weights))))),
    _normalVectors(std::move(normalVectors)),
    _surfaceLabel(surfaceLabel) {
  printInfo();
}

void MolecularSurfaceController::printInfo() {
  OutputControl::nOut << " --------------------------- " << std::endl;
  OutputControl::nOut << " Molecular Surface Information:" << std::endl;
  OutputControl::nOut << "  Surface Label:   " << _surfaceLabel << std::endl;
  OutputControl::nOut << "  N cavity points: " << _gridController->getNGridPoints() << std::endl;
  OutputControl::nOut << "  Surface area:    " << _gridController->getWeights().sum() << std::endl;
  OutputControl::nOut << " --------------------------- " << std::endl;
}

std::shared_ptr<GridController> MolecularSurfaceController::getGridController() {
  return _gridController;
}

const Eigen::Matrix3Xd& MolecularSurfaceController::getGridPoints() {
  return _gridController->getGridPoints();
}

const Eigen::VectorXd& MolecularSurfaceController::getWeights() {
  return _gridController->getWeights();
}

const Eigen::Matrix3Xd& MolecularSurfaceController::getNormalVectors() {
  return *_normalVectors;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixS() {
  if (!_S) {
    const double factor = _k * sqrt(4.0 * M_PI);
    unsigned int gridSize = _gridController->getNGridPoints();
    _S = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    Eigen::MatrixXd& s_inv = *_S;
    const Eigen::Matrix3Xd& coordinates = _gridController->getGridPoints();
    const Eigen::VectorXd& areas = _gridController->getWeights();
    for (unsigned int i = 0; i < gridSize; ++i) {
      s_inv(i, i) = factor / sqrt(areas(i));
      for (unsigned int j = 0; j < i; ++j) {
        double element = 1.0 / (coordinates.col(i) - coordinates.col(j)).norm();
        s_inv(i, j) = element;
        s_inv(j, i) = element;
      } // for j
    }   // for i
  }
  return *_S;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixSinv() {
  if (!_invS) {
    this->getMatrixS();
    _invS = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(_S->cols(), _S->cols()));
    *_invS = pseudoInvers_Sym(*_S);
  }
  return *_invS;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixD() {
  if (!_D) {
    unsigned int gridSize = _gridController->getNGridPoints();
    _D = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    Eigen::MatrixXd& d = *_D;
    const Eigen::VectorXd& a = _gridController->getWeights();
    const Eigen::Matrix3Xd& s = _gridController->getGridPoints();
    const Eigen::Matrix3Xd& n = *_normalVectors;
    //    const EigenVectorXd sphereRadii = *_sphereRadii;
    for (unsigned int i = 0; i < gridSize; ++i) {
      Eigen::Matrix3Xd coordDiff = s.colwise() - s.col(i);
      Eigen::VectorXd coordDiffNorm = coordDiff.colwise().norm().transpose();
      coordDiffNorm(i) = 1; // Avoid dividing by 0.
      Eigen::VectorXd d_i = (coordDiff.array() * n.array()).colwise().sum().matrix().transpose();
      d_i.array() /= (coordDiffNorm.array() * coordDiffNorm.array() * coordDiffNorm.array());
      // Overwrite diagonal entry.
      d_i(i) = -(2.0 * M_PI + (d_i.array() * a.array()).sum()) / a(i);
      d.col(i) = d_i;
    }
  }
  return *_D;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixA() {
  if (!_A) {
    unsigned int gridSize = _gridController->getNGridPoints();
    const Eigen::VectorXd& a = _gridController->getWeights();
    _A = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    _A->diagonal() = a;
  }
  return *_A;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixAinv() {
  if (!_Ainv) {
    unsigned int gridSize = _gridController->getNGridPoints();
    const Eigen::VectorXd& a = _gridController->getWeights();
    _Ainv = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    _Ainv->diagonal().array() = Eigen::VectorXd::Constant(gridSize, 1.0).array() / a.array();
  }
  return *_Ainv;
}

} /* namespace Serenity */
