/**
 * @file CoulombIntegralsOnGridController.cpp
 *
 * @author Moritz Bensberg
 * @date Aug 11, 2020
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
#include "data/grid/CoulombIntegralsOnGridController.h"
namespace Serenity {

CoulombIntegralsOnGridController::CoulombIntegralsOnGridController(std::shared_ptr<BasisController> basisController,
                                                                   std::shared_ptr<GridController> gridController,
                                                                   std::string fBaseName, unsigned int cacheSize,
                                                                   std::shared_ptr<Eigen::Matrix3Xd> normalVectors)
  : _basisController(basisController),
    _gridController(gridController),
    _intFileBaseName(fBaseName + ".elecPotInts.h5"),
    _nSetsCached(cacheSize),
    _normalVectors(normalVectors) {
  _gridController->addSensitiveObject(this->_self);
}

void CoulombIntegralsOnGridController::integralsToHDF5(std::vector<std::shared_ptr<Eigen::MatrixXd>> interals,
                                                       unsigned int blockIndex, HDF5::H5File& file) {
  std::string blockName = std::to_string(blockIndex);
  if (interals.size() > 0 && not _sizesInitialized) {
    _nRows = interals[0]->rows();
    _nCols = interals[0]->cols();
    _sizesInitialized = true;
  }
  const unsigned int nElements = _nRows * _nCols;
  const unsigned int blockSize = _blockSizes[blockIndex];
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ints(blockSize, nElements);
  for (unsigned int iSet = 0; iSet < blockSize; ++iSet) {
    Eigen::Map<Eigen::RowVectorXd> tmp(interals[iSet]->data(), nElements);
    ints.row(iSet) = tmp;
  } // for iSet
  HDF5::save(file, blockName, ints);
}

std::vector<std::shared_ptr<Eigen::MatrixXd>> CoulombIntegralsOnGridController::integralsFromHDF5(unsigned int blockIndex,
                                                                                                  HDF5::H5File& file) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ints;
  std::string blockName = std::to_string(blockIndex);
  HDF5::dataset_exists(file, blockName);
  HDF5::load(file, blockName, ints);
  const unsigned int blockSize = _blockSizes[blockIndex];
  std::vector<std::shared_ptr<Eigen::MatrixXd>> integralSets;
  for (unsigned int iSet = 0; iSet < blockSize; ++iSet) {
    auto newSet = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(_nRows, _nCols));
    Eigen::Map<Eigen::MatrixXd> tmp(ints.row(iSet).data(), _nRows, _nCols);
    *newSet = tmp;
    integralSets.push_back(newSet);
  }
  return integralSets;
}

Eigen::MatrixXd
CoulombIntegralsOnGridController::calculatePotentialIntegrals(std::vector<std::pair<double, std::array<double, 3>>> point) {
  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::nuclear, 0, 2, point, 0.0);
  Eigen::MatrixXd results = libint.compute1eInts(LIBINT_OPERATOR::nuclear, this->getBasisController(), point);
  libint.finalize(LIBINT_OPERATOR::nuclear, 0, 2);
  return results;
}

void CoulombIntegralsOnGridController::calculateElectricFieldIntegralsNumerically(
    Eigen::MatrixXd& result, std::vector<std::pair<double, std::array<double, 3>>> point, const Eigen::Vector3d& normalVector) {
  const double delta = 1e-6;
  const double factor = -1.0 / (delta);
  auto point_plus = point;
  auto point_minus = point;
  point_plus[0].second[0] += delta * normalVector(0);
  point_plus[0].second[1] += delta * normalVector(1);
  point_plus[0].second[2] += delta * normalVector(2);

  const Eigen::MatrixXd integralsPlus = this->calculatePotentialIntegrals(point_plus);
  const unsigned int nBFs = result.rows();
  result.rightCols(nBFs) = factor * (integralsPlus - result.leftCols(nBFs));
}

std::vector<std::shared_ptr<Eigen::MatrixXd>>
CoulombIntegralsOnGridController::calculateIntegralSet(const Eigen::Matrix3Xd& coordinates,
                                                       const Eigen::Matrix3Xd& normalVectors) {
  const unsigned int storeElecField = (_normalVectors) ? 2 : 1;
  const unsigned int nSets = coordinates.cols();
  std::vector<std::shared_ptr<Eigen::MatrixXd>> integrals(nSets, nullptr);

  const unsigned int nBFs = this->getBasisController()->getNBasisFunctions();
  Eigen::setNbThreads(1);
#pragma omp for schedule(dynamic)
  for (unsigned int iSet = 0; iSet < nSets; iSet++) {
    std::vector<std::pair<double, std::array<double, 3>>> point = {
        {-1.0, {{coordinates(0, iSet), coordinates(1, iSet), coordinates(2, iSet)}}}};
    integrals[iSet] = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(nBFs, storeElecField * nBFs));
    auto& results = *integrals[iSet];
    results.leftCols(nBFs) = calculatePotentialIntegrals(point);
    if (_normalVectors) {
      calculateElectricFieldIntegralsNumerically(results, point, normalVectors.col(iSet).eval());
    }
  } // gridPoint
  Eigen::setNbThreads(0);
  return integrals;
}

std::shared_ptr<BasisController> CoulombIntegralsOnGridController::getBasisController() {
  return _basisController;
}

void CoulombIntegralsOnGridController::initializeIntegrals() {
  const Eigen::Matrix3Xd& coordinates = _gridController->getGridPoints();
  const unsigned int nTotalSets = coordinates.cols();
  const unsigned int nBlocks = (unsigned int)ceil((double)nTotalSets / _nSetsCached);
  _cache = std::vector<std::shared_ptr<Eigen::MatrixXd>>(nTotalSets, nullptr);
  unsigned int start = 0;
  HDF5::H5File file(_intFileBaseName.c_str(), H5F_ACC_TRUNC);
  Libint::getInstance().keepEngines(LIBINT_OPERATOR::nuclear, 0, 2);
  if (_normalVectors)
    Libint::getInstance().keepEngines(LIBINT_OPERATOR::nuclear, 1, 2);
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned int n = _nSetsCached;
    if (iBlock == 0)
      n = nTotalSets - (nBlocks - 1) * _nSetsCached;
    _blockSizes.push_back(n);
    unsigned int blockEnd = start + n;
    std::vector<std::shared_ptr<Eigen::MatrixXd>> sets;
    sets = calculateIntegralSet(coordinates.block(0, start, 3, n),
                                (_normalVectors) ? _normalVectors->block(0, start, 3, n) : Eigen::Matrix3Xd(3, 1));

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

void CoulombIntegralsOnGridController::cleanUp() {
  std::string name = _intFileBaseName;
  std::remove(name.c_str());
  std::vector<unsigned int>().swap(_blockSizes);
  std::vector<std::shared_ptr<Eigen::MatrixXd>>().swap(_cache);
}

std::shared_ptr<GridController> CoulombIntegralsOnGridController::getGridController() {
  return _gridController;
}

void CoulombIntegralsOnGridController::notify() {
  cleanUp();
}

} /* namespace Serenity */
