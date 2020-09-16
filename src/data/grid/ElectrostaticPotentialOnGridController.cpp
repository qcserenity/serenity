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
#include "basis/BasisController.h"                      //getNBasisFunctions()
#include "data/grid/CoulombPotentialOnGridCalculator.h" //calculateElectronNuclei(), calculateElectronElectron()
#include "data/matrices/DensityMatrixController.h"      //getDensityMatrix()
#include "geometry/Geometry.h"                          //Nuclear potential
#include "integrals/wrappers/Libint.h"                  //compute()
#include "io/HDF5.h"                                    //Caching.
#include "misc/Timing.h"                                //Timings.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ElectrostaticPotentialOnGridController<SCFMode>::ElectrostaticPotentialOnGridController(
    std::shared_ptr<GridController> gridController, std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController,
    std::shared_ptr<Geometry> geometry, std::string fBaseName, unsigned int cacheSize)
  : _gridController(gridController),
    _densityMatrixController(densityMatrixController),
    _geometry(geometry),
    _fBaseName(fBaseName),
    _nSetsCached(cacheSize) {
  if (_densityMatrixController)
    _densityMatrixController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  _intFileBaseName = _fBaseName + ".elecPotInts.h5";
}

template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::cleanUpDisk() {
  std::string name = _intFileBaseName;
  std::remove(name.c_str());
  name = _fBaseName + ".elecPotGrid.h5";
  std::remove(name.c_str());
  _blockSizes.clear();
  _cache.clear();
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
FockMatrix<RESTRICTED>
ElectrostaticPotentialOnGridController<SCFMode>::integrateFockMatrix(const GridPotential<RESTRICTED>& charges) {
  assert(_densityMatrixController);
  assert(_electrostaticPotential->getGridController() == charges.getGridController());
  auto basisController = _densityMatrixController->getDensityMatrix().getBasisController();
  FockMatrix<RESTRICTED> f(basisController);
  f.setZero();
  const Eigen::Matrix3Xd& coordinates = _electrostaticPotential->getGridController()->getGridPoints();
  const unsigned int nTotalSets = coordinates.cols();
  if (nTotalSets < _nSetsCached)
    _nSetsCached = nTotalSets;
  if (_cache.size() == 0)
    initializeIntegrals(basisController, coordinates);
  const unsigned int nBlocks = _blockSizes.size();
  unsigned int blockEnd = nTotalSets;
  HDF5::Filepath name(_intFileBaseName);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  for (int iBlock = nBlocks - 1; iBlock >= 0; --iBlock) {
    unsigned int n = _blockSizes[iBlock];
    unsigned int start = blockEnd - n;
    if (!_cache[start]) {
      auto ints = integralsFromHDF5(iBlock, file);
      for (unsigned int i = start; i < blockEnd; ++i)
        _cache[i] = ints[i - start];
    }
    for (unsigned int i = start; i < blockEnd; ++i) {
      f -= charges[i] * (*_cache[i]);
      if (nBlocks > 1)
        _cache[i] = nullptr;
    }
    blockEnd = start;
  }
  file.close();
  return f;
}

template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::update() {
  const Eigen::Matrix3Xd& coordinates = _electrostaticPotential->getGridController()->getGridPoints();
  if (_densityMatrixController) {
    const unsigned int nTotalSets = coordinates.cols();
    if (nTotalSets < _nSetsCached)
      _nSetsCached = nTotalSets;
    if (_cache.size() == 0)
      initializeIntegrals(_densityMatrixController->getDensityMatrix().getBasisController(), coordinates);
    const unsigned int nBlocks = _blockSizes.size();
    unsigned int blockEnd = nTotalSets;
    auto& potential = *_electrostaticPotential;
    const auto totalDensityMatrix = _densityMatrixController->getDensityMatrix().total();
    HDF5::Filepath name(_intFileBaseName);
    HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    for (int iBlock = nBlocks - 1; iBlock >= 0; --iBlock) {
      unsigned int n = _blockSizes[iBlock];
      unsigned int start = blockEnd - n;
      if (!_cache[start]) {
        auto ints = integralsFromHDF5(iBlock, file);
        for (unsigned int i = start; i < blockEnd; ++i)
          _cache[i] = ints[i - start];
      }
      for (unsigned int i = start; i < blockEnd; ++i) {
        potential[i] = (_cache[i]->array() * totalDensityMatrix.array()).sum();
        if (nBlocks > 1)
          _cache[i] = nullptr;
      }
      blockEnd = start;
    }
    file.close();
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
  _electrostaticPotential = std::make_shared<GridPotential<RESTRICTED>>(_gridController);
  if (_diskMode && _diskUpToDate) {
    fromHDF5();
  }
  else {
    update();
  }
  Timings::timeTaken(" Tech. -             Elec. Pot.");
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<Eigen::MatrixXd>>
ElectrostaticPotentialOnGridController<SCFMode>::calculateIntegralSet(std::shared_ptr<BasisController> basisController,
                                                                      Eigen::Matrix3Xd coordinates) {
  const unsigned int nSets = coordinates.cols();
  std::vector<std::shared_ptr<Eigen::MatrixXd>> integrals(nSets, nullptr);

  auto& libint = Libint::getInstance();
  Eigen::setNbThreads(1);
#pragma omp for schedule(dynamic)
  for (unsigned int iSet = 0; iSet < nSets; iSet++) {
    std::vector<std::pair<double, std::array<double, 3>>> point = {
        {-1.0, {{coordinates(0, iSet), coordinates(1, iSet), coordinates(2, iSet)}}}};
    integrals[iSet] =
        std::make_shared<Eigen::MatrixXd>(libint.compute1eInts(libint2::Operator::nuclear, basisController, point));
  } // gridPoint
  Eigen::setNbThreads(0);
  return integrals;
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::initializeIntegrals(std::shared_ptr<BasisController> basisController,
                                                                          const Eigen::Matrix3Xd& coordinates) {
  const unsigned int nTotalSets = coordinates.cols();
  const unsigned int nBlocks = (unsigned int)ceil((double)nTotalSets / _nSetsCached);
  _cache = std::vector<std::shared_ptr<Eigen::MatrixXd>>(nTotalSets, nullptr);
  unsigned int start = 0;
  HDF5::H5File file(_intFileBaseName.c_str(), H5F_ACC_TRUNC);
  Libint::getInstance().keepEngines(libint2::Operator::nuclear, 0, 2);
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned int n = _nSetsCached;
    if (iBlock == 0)
      n = nTotalSets - (nBlocks - 1) * _nSetsCached;
    _blockSizes.push_back(n);
    unsigned int blockEnd = start + n;
    auto sets = calculateIntegralSet(basisController, coordinates.block(0, start, 3, n));
    if (nBlocks > 0)
      integralsToHDF5(sets, iBlock, file);
    // Keep last block in memory.
    if (iBlock == nBlocks - 1) {
      for (unsigned int i = start; i < blockEnd; ++i)
        _cache[i] = sets[i - start];
    }
    start += n;
  } // for iBlock
  file.close();
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
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::dataset_exists(file, "electrostaticPotential");
  HDF5::load(file, "electrostaticPotential", *_electrostaticPotential);
  file.close();
}
template<Options::SCF_MODES SCFMode>
void ElectrostaticPotentialOnGridController<SCFMode>::integralsToHDF5(std::vector<std::shared_ptr<Eigen::MatrixXd>> interals,
                                                                      unsigned int blockIndex, HDF5::H5File& file) {
  std::string blockName = std::to_string(blockIndex);
  const unsigned int nRows = _densityMatrixController->getDensityMatrix().getBasisController()->getNBasisFunctions();

  const unsigned int nElements = nRows * nRows;
  const unsigned int blockSize = _blockSizes[blockIndex];
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ints(blockSize, nElements);
  for (unsigned int iSet = 0; iSet < blockSize; ++iSet) {
    Eigen::Map<Eigen::RowVectorXd> tmp(interals[iSet]->data(), nElements);
    ints.row(iSet) = tmp;
  } // for iSet
  HDF5::save(file, blockName, ints);
  //  file.close();
}
template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<Eigen::MatrixXd>>
ElectrostaticPotentialOnGridController<SCFMode>::integralsFromHDF5(unsigned int blockIndex, HDF5::H5File& file) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ints; // = std::make_shared<Eigen::MatrixXd>();
  std::string blockName = std::to_string(blockIndex);
  HDF5::dataset_exists(file, blockName);
  HDF5::load(file, blockName, ints);
  const unsigned int blockSize = _blockSizes[blockIndex];
  std::vector<std::shared_ptr<Eigen::MatrixXd>> integralSets;
  const unsigned int nRows = _densityMatrixController->getDensityMatrix().getBasisController()->getNBasisFunctions();
  for (unsigned int iSet = 0; iSet < blockSize; ++iSet) {
    auto newSet = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(nRows, nRows));
    Eigen::Map<Eigen::MatrixXd> tmp(ints.row(iSet).data(), nRows, nRows);
    *newSet = tmp;
    integralSets.push_back(newSet);
  }
  return integralSets;
}

template class ElectrostaticPotentialOnGridController<Options::SCF_MODES::RESTRICTED>;
template class ElectrostaticPotentialOnGridController<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
