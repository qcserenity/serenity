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
/* Include Class Header*/
#include "geometry/MolecularSurfaceController.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"                  //addSensitiveObject()
#include "geometry/MolecularSurface.h"          //Underlying surface
#include "geometry/MolecularSurfaceFactory.h"   //Surface construction.
#include "geometry/Sphere.h"                    //Cavity energy.
#include "io/FormattedOutputStream.h"           //printInfo()
#include "math/linearAlgebra/MatrixFunctions.h" //Pseudo inversion.
#include "misc/WarningTracker.h"                //Warnings.
#include "parameters/Constants.h"               //Gas constant.
#include "settings/PCMSettings.h"               //Get the number density.
#include "solvation/Solvents.h"                 //The tabulated solvent data.
/* Include Std and External Headers */
#include <cmath> //M_PI and log.

namespace Serenity {

MolecularSurfaceController::MolecularSurfaceController(std::shared_ptr<Geometry> geometry, const PCMSettings& pcmSettings)
  : _geometry(geometry), _pcmSettings(pcmSettings) {
  _geometry->addSensitiveObject(this->_self);
}

const MolecularSurface& MolecularSurfaceController::getMolecularSurface() {
  /*
   * This cast looks rather ugly and I have not added any error-handling for
   * a potential failure. However, the _grid which is accessed and casted can/
   * should only be constructed in this class as a MolecularSurface. Thus, this
   * cast should always work. I have implemented it like this, since the
   * MolecularSurfaceController is a GridController which makes it
   * necessary to hold _grid as a unique_ptr.
   */
  return dynamic_cast<MolecularSurface&>(*this->_grid);
}

void MolecularSurfaceController::printInfo() {
  if (!this->_grid)
    buildSurface();
  auto& surface = getMolecularSurface();
  OutputControl::nOut << " -------------------------------- " << std::endl;
  OutputControl::nOut << " Molecular Surface Information:" << std::endl;
  OutputControl::nOut << "  Surface Label:   " << surface.getLabel() << std::endl;
  OutputControl::nOut << "  N cavity points: " << this->getNGridPoints() << std::endl;
  OutputControl::nOut << "  Surface area:    " << this->getWeights().sum() << std::endl;
  OutputControl::nOut << " -------------------------------- " << std::endl;
}

const Eigen::Matrix3Xd& MolecularSurfaceController::getNormalVectors() {
  if (!this->_grid)
    buildSurface();
  return this->getMolecularSurface().getNormalVectors();
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixS() {
  if (!this->_grid)
    buildSurface();
  if (!_S) {
    double minDist = 10;
    const double factor = _k * sqrt(4.0 * M_PI);
    unsigned int gridSize = this->getNGridPoints();
    _S = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    Eigen::MatrixXd& s_inv = *_S;
    const Eigen::Matrix3Xd& coordinates = this->getGridPoints();
    const Eigen::VectorXd& areas = this->getWeights();
    for (unsigned int i = 0; i < gridSize; ++i) {
      s_inv(i, i) = factor / sqrt(areas(i));
      for (unsigned int j = 0; j < i; ++j) {
        double element = 1.0 / (coordinates.col(i) - coordinates.col(j)).norm();
        minDist = std::min(1.0 / element, minDist);
        s_inv(i, j) = element;
        s_inv(j, i) = element;
      } // for j
    }   // for i
    OutputControl::vOut << "Mol.Surface: Min. inter cavity-point distance: " << minDist << std::endl;
  }
  return *_S;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixSinv() {
  if (!this->_grid)
    buildSurface();
  if (!_invS) {
    this->getMatrixS();
    _invS = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(_S->cols(), _S->cols()));
    *_invS = pseudoInvers_Sym(*_S);
    double maxValue = _invS->array().abs().maxCoeff();
    OutputControl::vOut << "Mol.Surface: Abs max. value of S^(-1): " << maxValue << std::endl;
    if (maxValue > 10) {
      WarningTracker::printWarning(
          "WARNING: Large value in S^(-1) detected. This may lead to unstable gradient/SCF calculations.\n"
          "         Consider increasing <minDistance> in the PCM input-block.",
          true);
    }
  }
  return *_invS;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixD() {
  if (!this->_grid)
    buildSurface();
  if (!_D) {
    unsigned int gridSize = this->getNGridPoints();
    _D = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    Eigen::MatrixXd& d = *_D;
    const Eigen::VectorXd& a = this->getWeights();
    const Eigen::Matrix3Xd& s = this->getGridPoints();
    const Eigen::Matrix3Xd& n = this->getNormalVectors();
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
  if (!this->_grid)
    buildSurface();
  if (!_A) {
    unsigned int gridSize = this->getNGridPoints();
    const Eigen::VectorXd& a = this->getWeights();
    _A = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    _A->diagonal() = a;
  }
  return *_A;
}

const Eigen::MatrixXd& MolecularSurfaceController::getMatrixAinv() {
  if (!this->_grid)
    buildSurface();
  if (!_Ainv) {
    unsigned int gridSize = this->getNGridPoints();
    const Eigen::VectorXd& a = this->getWeights();
    _Ainv = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(gridSize, gridSize));
    _Ainv->diagonal().array() = Eigen::VectorXd::Constant(gridSize, 1.0).array() / a.array();
  }
  return *_Ainv;
}

void MolecularSurfaceController::buildSurface() {
  _grid = MolecularSurfaceFactory::produce(_geometry, _pcmSettings);
  printInfo();
}

const std::vector<unsigned int>& MolecularSurfaceController::getPointToSphereMapping() {
  if (!this->_grid)
    buildSurface();
  return getMolecularSurface().pointWiseSphereIndices();
}

const std::vector<std::pair<unsigned int, unsigned int>>& MolecularSurfaceController::getSphereToPointMapping() {
  if (!this->_grid)
    buildSurface();
  return getMolecularSurface().getGridIndicesOfAtoms();
}

void MolecularSurfaceController::notify() {
  if (this->isLoaded()) // if pcm loaded, give warning
    WarningTracker::printWarning("Warning: A molecular surface for solvation was loaded, but now the molecular surface "
                                 "is being reset! The loaded surface may not be used as intended and may now change!",
                                 true);
  _grid = nullptr;
  _S = nullptr;
  _invS = nullptr;
  _D = nullptr;
  _A = nullptr;
  _Ainv = nullptr;
  this->notifyObjects();
}

const Eigen::Matrix3Xd& MolecularSurfaceController::getGridPoints() {
  if (!this->_grid)
    buildSurface();
  return this->GridController::getGridPoints();
}

const Eigen::VectorXd& MolecularSurfaceController::getWeights() {
  if (!this->_grid)
    buildSurface();
  return this->GridController::getWeights();
}

unsigned int MolecularSurfaceController::getNGridPoints() {
  if (!this->_grid)
    buildSurface();
  return this->GridController::getNGridPoints();
}

std::shared_ptr<Geometry> MolecularSurfaceController::getGeometry() {
  return _geometry;
}

double MolecularSurfaceController::getCavityEnergy() {
  if (!_cavityEnergy)
    calculateCavityEnergy();
  return *_cavityEnergy;
}

void MolecularSurfaceController::setSurface(std::unique_ptr<MolecularSurface>&& surface) {
  this->_grid = std::move(surface);
}

bool MolecularSurfaceController::isLoaded() {
  return this->_pcmSettings.loadedPCM;
}

std::string MolecularSurfaceController::getChargesPath() {
  return this->_pcmSettings.cavityPath;
}

void MolecularSurfaceController::calculateCavityEnergy() {
  if (!this->_grid) {
    OutputControl::nOut << " Auxiliary Cavity Construction" << std::endl;
    buildSurface();
  }
  const auto& sphereToPointMap = this->getMolecularSurface().getGridIndicesOfAtoms();
  const double nSpheres = sphereToPointMap.size();
  const Eigen::VectorXd& weights = this->getWeights();
  const double r_s = this->getCavityFormationProbeRadius();
  const auto& spheres = this->getMolecularSurface().getSpheres();
  const double numberDensity = this->getNumberDensity();
  _cavityEnergy = std::make_shared<double>();
  *_cavityEnergy = 0.0;
  if (nSpheres != spheres.size())
    throw SerenityError(
        (std::string) "ERROR: The number of spheres used in the cavity representation is inconsistent!\n" +
        "       This may be caused by using a not fully supported surface type. Please use\n" +
        "       The Delley-type surface.");
  OutputControl::nOut << " ------------------------------------------------------ " << std::endl;
  OutputControl::nOut << " Cavity formation energy:" << std::endl;
  OutputControl::nOut << "  Number density (1/bohr^3):   " << numberDensity << std::endl;
  OutputControl::nOut << "  Probe radius:                " << r_s << std::endl;
  OutputControl::nOut << "  Temperature:                 " << _pcmSettings.temperature << std::endl;
  for (unsigned int iSphere = 0; iSphere < nSpheres; ++iSphere) {
    const double a = calculateExposedSurface_sphere(sphereToPointMap[iSphere], weights);
    const double r_m = spheres[iSphere].getRadius();
    const double g_cav_sphere = calculateG_cav_sphere(r_s, numberDensity, r_m, _pcmSettings.temperature);
    *_cavityEnergy += a / (4.0 * M_PI * r_m * r_m) * g_cav_sphere;
  };
  OutputControl::nOut << "  Cavity Formation Energy:     " << *_cavityEnergy << std::endl;
  OutputControl::nOut << " ------------------------------------------------------ " << std::endl;
}

double MolecularSurfaceController::calculateG_cav_sphere(double r_s, double rho, double r_m, double t) {
  const double r_s2 = 2.0 * r_s;
  const double y = M_PI / 6.0 * r_s2 * r_s2 * r_s2 * rho;
  const double rm_rs = r_m / r_s;
  const double rt = UNIVERSALGAS_CONSTANT * t / (1000 * HARTREE_TO_KJ_PER_MOL);
  const double y_1my = y / (1.0 - y);
  double g_cav = -log(1.0 - y) + 3.0 * y_1my * rm_rs + (3.0 * y_1my + 9.0 / 2.0 * y_1my * y_1my) * rm_rs * rm_rs;
  g_cav *= rt;
  return g_cav;
}

double MolecularSurfaceController::calculateExposedSurface_sphere(const std::pair<unsigned int, unsigned int> indices,
                                                                  const Eigen::VectorXd& weights) {
  int nPoints = indices.second - indices.first;
  if (nPoints <= 0)
    return 0.0;
  return weights.segment(indices.first, nPoints).sum();
}

double MolecularSurfaceController::getNumberDensity() {
  if (_pcmSettings.solvent == Options::PCM_SOLVENTS::EXPLICIT)
    return _pcmSettings.numberDensity;
  return Solvents::getNumberDensity(_pcmSettings.solvent);
}

double MolecularSurfaceController::getCavityFormationProbeRadius() {
  if (_pcmSettings.solvent == Options::PCM_SOLVENTS::EXPLICIT)
    return _pcmSettings.cavityProbeRadius;
  return Solvents::getCavityFormationProbeRadius(_pcmSettings.solvent);
}

} /* namespace Serenity */
