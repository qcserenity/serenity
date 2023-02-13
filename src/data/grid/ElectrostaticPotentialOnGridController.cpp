/**
 * @file ElectrostaticPotentialOnGridController.cpp
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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
#include "data/grid/ElectrostaticPotentialOnGridController.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"          //Basis function to atom mapping.
#include "basis/BasisController.h"                      //getNBasisFunctions()
#include "data/grid/CoulombPotentialOnGridCalculator.h" //calculateElectronNuclei(), calculateElectronElectron()
#include "data/matrices/DensityMatrixController.h"      //getDensityMatrix()
#include "geometry/Geometry.h"                          //Nuclear potential
#include "integrals/wrappers/Libint.h"                  //compute()
#include "io/HDF5.h"                                    //Caching.
#include "misc/SerenityError.h"                         //Throw errors.
#include "misc/Timing.h"                                //Timings.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ElectrostaticPotentialOnGridController<SCFMode>::ElectrostaticPotentialOnGridController(
    std::shared_ptr<GridController> gridController,
    std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController, std::shared_ptr<Geometry> geometry,
    std::string fBaseName, unsigned int cacheSize, std::shared_ptr<Eigen::Matrix3Xd> normalVectors,
    std::shared_ptr<AtomCenteredBasisController> atomCenteredBasisController)
  : _coulombIntegralsOnDiskController(std::make_shared<CoulombIntegralsOnGridController>(
        densityMatrixController->getDensityMatrix().getBasisController(), gridController, fBaseName, cacheSize, normalVectors)),
    _densityMatrixController(densityMatrixController),
    _geometry(geometry),
    _fBaseName(fBaseName),
    _atomCenteredBasisController(atomCenteredBasisController) {
  _densityMatrixController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  _coulombIntegralsOnDiskController->getGridController()->addSensitiveObject(ObjectSensitiveClass<Grid>::_self);
}

template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::cleanUpDisk() {
  _coulombIntegralsOnDiskController->cleanUp();
  std::string name = _fBaseName + ".elecPotGrid.h5";
  std::remove(name.c_str());
  _diskUpToDate = false;
}

template<Options::SCF_MODES SCFMode>
ElectrostaticPotentialOnGridController<SCFMode>::~ElectrostaticPotentialOnGridController() {
  cleanUpDisk();
}

template<Options::SCF_MODES SCFMode>
const GridPotential<RESTRICTED>& ElectrostaticPotentialOnGridController<SCFMode>::getPotential() {
  if (!_electrostaticPotential)
    getData();
  return *_electrostaticPotential;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<BasisController> ElectrostaticPotentialOnGridController<SCFMode>::getBasisController() {
  return _coulombIntegralsOnDiskController->getBasisController();
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<GridController> ElectrostaticPotentialOnGridController<SCFMode>::getGridController() {
  return _coulombIntegralsOnDiskController->getGridController();
}

template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::update() {
  unsigned int nColsToContract = getBasisController()->getNBasisFunctions();
  if (_densityMatrixController) {
    auto potential = _electrostaticPotential;
    const auto totalDensityMatrix = _densityMatrixController->getDensityMatrix().total();
    auto func = [potential, nColsToContract, &totalDensityMatrix](std::shared_ptr<Eigen::MatrixXd> ints, unsigned int i,
                                                                  unsigned int /*tid*/) {
      potential->operator[](i) = (ints->leftCols(nColsToContract).array() * totalDensityMatrix.array()).sum();
    };
    _coulombIntegralsOnDiskController->iterateGridPoints(func);
  }
  if (_geometry)
    CoulombPotentialOnGridCalculator::calculateElectronNuclei(*_electrostaticPotential, _geometry->getAtoms());
  // switch sign, since the CoulombPotentialOnGridCalculator has the sign inverted.
  *_electrostaticPotential *= -1.0;
  if (_diskMode)
    toHDF5();
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::getData() {
  Timings::takeTime(" Tech. -             Elec. Pot.");
  _electrostaticPotential = std::make_shared<GridPotential<RESTRICTED>>(this->getGridController());
  if (_diskMode && _diskUpToDate) {
    fromHDF5();
  }
  else {
    update();
  }
  Timings::timeTaken(" Tech. -             Elec. Pot.");
}

template<Options::SCF_MODES SCFMode>
FockMatrix<RESTRICTED>
ElectrostaticPotentialOnGridController<SCFMode>::integrateFockMatrix(const GridPotential<RESTRICTED>& charges) {
  assert(this->getGridController() == charges.getGridController());
  FockMatrix<RESTRICTED> f(this->getBasisController());
  f.setZero();
  unsigned int nColsToContract = this->getBasisController()->getNBasisFunctions();
  unsigned int nThreads = omp_get_max_threads();
  std::vector<Eigen::MatrixXd> threadWiseMatrices(nThreads, Eigen::MatrixXd::Zero(nColsToContract, nColsToContract));
  auto func = [&](std::shared_ptr<Eigen::MatrixXd> ints, unsigned int i, unsigned int threadId) {
    threadWiseMatrices[threadId] -= charges[i] * ints->leftCols(nColsToContract);
  };
  _coulombIntegralsOnDiskController->iterateGridPoints(func);
  for (const auto& mat : threadWiseMatrices)
    f += mat;
  return f;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd
ElectrostaticPotentialOnGridController<SCFMode>::calculateGradient(const std::vector<unsigned int>& atomIndicesOfPoints,
                                                                   const GridPotential<RESTRICTED>& charges) {
  // TODO: Introduce atom mapping. The geometry for the potential does not have to be the same as for the active
  // system/basis.
  const unsigned int nAtoms = _geometry->getNAtoms();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(nAtoms, 3);
  if (_densityMatrixController) {
    const auto totalDensityMatrix = _densityMatrixController->getDensityMatrix().total();
    unsigned int nThreads = omp_get_max_threads();
    std::vector<Eigen::MatrixXd> threadWiseMatrices(nThreads, Eigen::MatrixXd::Zero(nAtoms, 3));
    if (!_atomCenteredBasisController) {
      throw SerenityError("An atom centered basis controller is required to calculate gradients!");
    }
    const auto& atomIndicesOfBasis = _atomCenteredBasisController->getAtomIndicesOfBasis();
    auto func = [&](Eigen::MatrixXd& ints, unsigned int iPoint, unsigned int ii, unsigned int jj, unsigned int nI,
                    unsigned int nJ, unsigned int threadId) {
      Eigen::MatrixXd p_ij = totalDensityMatrix.block(ii, jj, nI, nJ);
      // Take care of the double occurrence of basis function combinations.
      // For identical shells (ii==jj) everything is already correct.
      const double perm = (ii == jj) ? 1.0 : 2.0;
      // d/dr chi_i op chi_j
      const double dx_i = (Eigen::Map<Eigen::MatrixXd>(ints.col(0).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dy_i = (Eigen::Map<Eigen::MatrixXd>(ints.col(1).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dz_i = (Eigen::Map<Eigen::MatrixXd>(ints.col(2).data(), nI, nJ).array() * p_ij.array()).sum();
      // chi_i op d/dr chi_j
      const double dx_j = (Eigen::Map<Eigen::MatrixXd>(ints.col(3).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dy_j = (Eigen::Map<Eigen::MatrixXd>(ints.col(4).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dz_j = (Eigen::Map<Eigen::MatrixXd>(ints.col(5).data(), nI, nJ).array() * p_ij.array()).sum();
      threadWiseMatrices[threadId](atomIndicesOfBasis[ii], 0) += perm * charges[iPoint] * dx_i;
      threadWiseMatrices[threadId](atomIndicesOfBasis[ii], 1) += perm * charges[iPoint] * dy_i;
      threadWiseMatrices[threadId](atomIndicesOfBasis[ii], 2) += perm * charges[iPoint] * dz_i;
      threadWiseMatrices[threadId](atomIndicesOfBasis[jj], 0) += perm * charges[iPoint] * dx_j;
      threadWiseMatrices[threadId](atomIndicesOfBasis[jj], 1) += perm * charges[iPoint] * dy_j;
      threadWiseMatrices[threadId](atomIndicesOfBasis[jj], 2) += perm * charges[iPoint] * dz_j;
      // chi_i chi_j d/dr op
      const double dx_opt = (Eigen::Map<Eigen::MatrixXd>(ints.col(6).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dy_opt = (Eigen::Map<Eigen::MatrixXd>(ints.col(7).data(), nI, nJ).array() * p_ij.array()).sum();
      const double dz_opt = (Eigen::Map<Eigen::MatrixXd>(ints.col(8).data(), nI, nJ).array() * p_ij.array()).sum();
      threadWiseMatrices[threadId](atomIndicesOfPoints[iPoint], 0) += perm * charges[iPoint] * dx_opt;
      threadWiseMatrices[threadId](atomIndicesOfPoints[iPoint], 1) += perm * charges[iPoint] * dy_opt;
      threadWiseMatrices[threadId](atomIndicesOfPoints[iPoint], 2) += perm * charges[iPoint] * dz_opt;
    };
    _coulombIntegralsOnDiskController->iterateGridPoints_gradient(func);
    for (const auto& mat : threadWiseMatrices)
      gradient -= mat;
  }
  return gradient;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<DensityMatrixController<SCFMode>> ElectrostaticPotentialOnGridController<SCFMode>::getDensityMatrixController() {
  return _densityMatrixController;
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::setDiskMode(bool newMode) {
  _diskMode = newMode;
  if (_diskMode) {
    toHDF5();
    _electrostaticPotential.reset();
  }
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::toHDF5() {
  if (!_electrostaticPotential)
    update();
  std::string name = _fBaseName + ".elecPotGrid.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "electrostaticPotential", *_electrostaticPotential);
  file.close();
  _diskUpToDate = true;
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::fromHDF5() {
  HDF5::Filepath name(_fBaseName + ".elecPotGrid.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "electrostaticPotential");
  HDF5::load(file, "electrostaticPotential", *_electrostaticPotential);
  file.close();
}

template class ElectrostaticPotentialOnGridController<Options::SCF_MODES::RESTRICTED>;
template class ElectrostaticPotentialOnGridController<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
