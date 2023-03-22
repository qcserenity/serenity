/**
 * @file   CDStorageController.cpp
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
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
#include "integrals/CDStorageController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "integrals/CDIntegralController.h"
#include "memory/MemoryManager.h"
/* Include Std and External Headers */
#include <H5Cpp.h>
#include <hdf5.h>
#if __linux__ || __unix__ || __unix
#include <malloc.h> //Free unused memory.
#endif

namespace Serenity {

CDStorageController::CDStorageController(std::string systemName, std::string label,
                                         std::shared_ptr<CDIntegralController> cdIntController)
  : _systemName(systemName),
    _diskMode(cdIntController->getDiskMode()),
    _cdIntegralController(cdIntController),
    _vectorPool(new std::vector<std::shared_ptr<Eigen::RowVectorXd>>),
    _label(label),
    _doubletype(H5Tcopy(H5T_NATIVE_DOUBLE)),
    _vectorDim(0),
    _lastBatchStart(0),
    _lastBatchSize(0),
    _vectorsAllocated(false),
    _vectorsUpToDate(false) {
  this->addMatrix();
}

CDStorageController::~CDStorageController() {
  try {
    _file->close();
  }
  catch (...) {
  }
  std::remove((_systemName + ".cd." + _label + ".h5").c_str());
  H5Tclose(_doubletype);
};

void CDStorageController::addMatrix() {
  int nelemts;   /* Dummy parameter in API, no longer used */
  size_t nslots; /* Number of slots in the hash table */
  size_t nbytes; /* Size of chunk cache in bytes */
  double w0;     /* Chunk preemption policy */
  H5::FileAccPropList fapl;
  fapl.getCache(nelemts, nslots, nbytes, w0);
  nbytes = 2 * 1024 * 1024;
  fapl.setCache(nelemts, nslots, nbytes, w0);

  std::string name = _systemName + ".cd." + _label + ".h5";
  try {
    H5::H5File file(name.c_str(), H5F_ACC_TRUNC);
    _file = std::make_shared<H5::H5File>(file);
    _file->openFile(name.c_str(), H5F_ACC_RDWR, fapl);
  }
  catch (H5::FileIException& file_exists_err) {
    throw std::runtime_error(
        "A HDF5 file for the system is already open or has not been closed correctly! The File is: " + _systemName +
        ".cd." + _label + ".h5");
  }
}

void CDStorageController::freeLastBatch() {
  if (*_diskMode) {
    for (unsigned int index = _lastBatchStart; index < _lastBatchStart + _lastBatchSize; index++) {
      (*_vectorPool)[index].reset();
    }
  }
}

void CDStorageController::storeVector(const int index, std::shared_ptr<Eigen::RowVectorXd> vector) {
  _vectorPool->push_back(nullptr);
  if (*_diskMode) {
    HDF5::save((*_file), std::to_string(index), (*vector));
    this->updateDataSets();
  }
  if (!(*_diskMode)) {
    (*_vectorPool)[index] = vector;
    auto memManager = MemoryManager::getInstance();
    if (memManager->getAvailableSystemMemory() < 4e+9) {
      auto cdIntegralController = _cdIntegralController.lock();
      // Turn DiskMode on, save all vectors to disk, free memory
      cdIntegralController->setDiskMode();
    }
  }
}

void CDStorageController::storeDiag(std::shared_ptr<Eigen::VectorXd> vector) {
  HDF5::save((*_file), "diagonal", (*vector));
  _vectorDim = vector->size();
}

std::shared_ptr<Eigen::VectorXd> CDStorageController::loadDiag() {
  Eigen::VectorXd tmp(_vectorDim);
  H5::DataSet dataSet = _file->openDataSet("diagonal");
  HDF5::dataset_exists((*_file), "diagonal");
  HDF5::loadDirectly(dataSet, tmp);
  return std::make_shared<Eigen::VectorXd>(tmp);
}

void CDStorageController::updateDataSets() {
  unsigned int counter = _vectorPool->size();
  int size = _dataSets.size();
  for (unsigned int i = size; i < counter; i++) {
    _dataSets.push_back(_file->openDataSet(std::to_string(i)));
  }
}

void CDStorageController::setVectorDimension(unsigned int dim) {
  _vectorDim = dim;
}

unsigned int CDStorageController::getNVectors() {
  return _vectorPool->size();
}

std::shared_ptr<Eigen::RowVectorXd> CDStorageController::loadVector(const int index) {
  if (!(*_vectorPool)[index]) {
    throw SerenityError("CDStorageController: Vector not in Memory. Use loadBatch() before loadVector()");
  }
  return (*_vectorPool)[index];
}

// load vectors in batch from disk
unsigned int CDStorageController::loadBatchDisk(unsigned int startIndex, unsigned int batchSize) {
  Timings::takeTime("Chol. - load Vectors from Disk");
  auto memManager = MemoryManager::getInstance();
  auto availableMem = memManager->getAvailableSystemMemory();
  long long maxBatchSize = std::floor((availableMem - 3e+9) / (_vectorDim * 8));
  if (batchSize > maxBatchSize)
    batchSize = maxBatchSize;
  if (batchSize + startIndex > this->getNVectors())
    batchSize = this->getNVectors() - startIndex;
  _lastBatchSize = batchSize;
  _lastBatchStart = startIndex;
  for (unsigned int index = _lastBatchStart; index < _lastBatchStart + _lastBatchSize; index++) {
    if (!(*_vectorPool)[index]) {
      Eigen::RowVectorXd tmp(_vectorDim);
      HDF5::loadDirectly(_dataSets[index], tmp);
      (*_vectorPool)[index] = std::make_shared<Eigen::RowVectorXd>(tmp);
    }
  }
  Timings::timeTaken("Chol. - load Vectors from Disk");
  return _lastBatchSize;
}

// return the number of all vectors kept in memory
unsigned int CDStorageController::loadBatchMem(unsigned int startIndex, unsigned int batchSize) {
  _lastBatchStart = startIndex;
  assert(_lastBatchStart == 0);
  if (batchSize + startIndex > _vectorPool->size())
    batchSize = _vectorPool->size() - startIndex;
  _lastBatchSize = batchSize;
  return batchSize;
}

unsigned int CDStorageController::loadBatch(unsigned int startIndex, unsigned int batchSize) {
  if (*_diskMode) {
    if (_lastBatchStart == startIndex) {
      if ((*_vectorPool)[startIndex])
        return _lastBatchSize;
    }
    this->freeLastBatch();
    return loadBatchDisk(startIndex, batchSize);
  }
  return loadBatchMem(startIndex, batchSize);
}

void CDStorageController::freeMemory() {
  if (_vectorPool->size() == 0)
    return;
  for (unsigned int i = 0; i < _vectorPool->size(); i++) {
    // Free memory used to hold Vectors in Memory
    (*_vectorPool)[i].reset();
  }
  this->updateDataSets();
}

void CDStorageController::setDiskMode() {
  if (_vectorPool->size() == 0)
    return;
  for (unsigned int i = 0; i < _vectorPool->size(); i++) {
    // Write Vectors from Memory to Disk
    if ((*_vectorPool)[i]) {
      HDF5::save((*_file), std::to_string(i), *((*_vectorPool)[i]));
    }
    // Free memory used to hold Vectors in Memory
    (*_vectorPool)[i].reset();
  }
  this->updateDataSets();
  *(_diskMode) = true;
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

void CDStorageController::flushFile() {
  _file->flush(H5F_SCOPE_GLOBAL);
}

void CDStorageController::allocateVectors(unsigned int nVec, unsigned int dim) {
  _vectorPool->clear();

  // leave 2 GB of memory to ensure stability
  double memBuffer = 2e+9;
  auto memManager = MemoryManager::getInstance();
  double availableMem = memManager->getAvailableSystemMemory();
  //  availableMem -= 4e+9;
  availableMem -= sizeof(Eigen::RowVectorXd) * nVec;
  availableMem -= sizeof(std::make_shared<Eigen::RowVectorXd>(dim)) * nVec;
  availableMem -= double(sizeof(double) * double(nVec) * double(dim));

  if (availableMem > memBuffer) {
    *_diskMode = false;
    // allocate in mem
    while (_vectorPool->size() < nVec) {
      _vectorPool->push_back(std::make_shared<Eigen::RowVectorXd>(dim));
      _vectorPool->back()->setZero();
    }
  }
  else {
    // allocate on disk
    *_diskMode = true;
    // create contiguous datasets. Should allocate disk space on creation
    // rank of the datasets
    unsigned int rank = 1;
    // dataset dimension on file
    hsize_t dimsf[rank];
    // Creat dataspace
    dimsf[0] = dim;
    H5::DataSpace dataspace(rank, dimsf, NULL);
    // DSet creation property list
    H5::DSetCreatPropList plist = H5::DSetCreatPropList::DEFAULT;
    plist.setAllocTime(H5D_ALLOC_TIME_EARLY);
    // Creat dataset
    for (unsigned int i = 0; i < nVec; i++) {
      std::string name = std::to_string(i);
      _file->createDataSet(name.c_str(), _doubletype, dataspace, plist);
      // add placeholder to the vector pool
      _vectorPool->push_back(nullptr);
    }
    // update the vector holding all datasets
    this->updateDataSets();
  }
}

void CDStorageController::storeSegment(unsigned int index, unsigned int firstIndex,
                                       std::shared_ptr<std::vector<double>> val, unsigned int stride) {
  if (*_diskMode) {
    hsize_t dim[1];
    dim[0] = val->size();
    hsize_t strides[1];
    strides[0] = stride;
    double data[dim[0]];
    std::move(val->begin(), val->end(), data);

    auto dspace = _dataSets[index].getSpace();
    // dataset dimension on file
    hsize_t offset[1];
    offset[0] = firstIndex;

    dspace.selectHyperslab(H5S_SELECT_SET, dim, offset, strides);

    auto memspace = _dataSets[index].getSpace();
    // dataset dimension on file
    hsize_t moffset[1];
    moffset[0] = 0;

    memspace.selectHyperslab(H5S_SELECT_SET, dim, moffset, strides);

    _dataSets[index].write(data, _doubletype, memspace, dspace);

    dspace.close();
    memspace.close();
  }
  else {
    // This should be adjustable by stride
    constexpr unsigned int cStride = 1;
    // Left part of the equation:
    // Map the vector elements to a new vector object. The offset(firstIndex) is simply added to the data pointer,
    // while the stride is found in Eigen::InnerStride.
    // Right part of the equation:
    // The data of a std::vector<double> is mapped to an Eigen::VectorXd.
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<cStride>>((*_vectorPool)[index]->data() + firstIndex, val->size()) =
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(val->data(), val->size());
  }
}

void CDStorageController::notify() {
  if (_vectorDim > 0 || _vectorPool->size() > 0) {
    // remove old file from disk
    try {
      _file->close();
    }
    catch (...) {
    }
    std::remove((_systemName + ".cd." + _label + ".h5").c_str());

    // reset member variables
    auto cdIntCont = _cdIntegralController.lock();
    _file = nullptr;
    _diskMode = cdIntCont->getDiskMode();
    _vectorPool->clear();
    _dataSets.clear();
    _vectorDim = 0;
    _lastBatchStart = 0;
    _lastBatchSize = 0;
    _vectorsAllocated = false;
    _vectorsUpToDate = false;

    // setup a new file
    addMatrix();
  }
}

void CDStorageController::addSensitiveBasis(std::shared_ptr<BasisController> basCont) {
  basCont->addSensitiveObject(this->_self);
}

void CDStorageController::setUpToDate() {
  _vectorsUpToDate = true;
}

bool CDStorageController::getUpToDate() {
  return _vectorsUpToDate;
}

} /* namespace Serenity */
