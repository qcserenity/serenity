/**
 * @file VectorOnDiskStorageController.cpp
 *
 * @date Apr 19, 2018
 * @author Moritz Bensberg
 *
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
#include "misc/VectorOnDiskStorageController.h"
/* Include Serenity Internal Headers */
#include "io/Filesystem.h"
#include "io/HDF5.h"

namespace Serenity {

VectorOnDiskStorageController::VectorOnDiskStorageController(double cacheLimit, std::string fileName)
  : _cacheLimit(cacheLimit), _nDoubles(0), _fileName(fileName) {
}

VectorOnDiskStorageController::VectorOnDiskStorageController(VectorOnDiskStorageController& old, std::string fileName) {
  _cacheLimit = old.getCacheLimit();
  _fileName = fileName;
  _nDoubles = 0;
  const auto labelList = old.getLabelList();
  for (const auto& label : labelList) {
    this->storeVectorSegment(std::make_shared<Eigen::VectorXd>(*old.getVectorSegment(label)), label);
  }
}

std::vector<std::string> VectorOnDiskStorageController::getLabelList() {
  std::vector<std::string> labelList;
  for (const auto& it : _vectorPool) {
    labelList.push_back(it.first);
  }
  return labelList;
}

void VectorOnDiskStorageController::storeVectorSegment(std::shared_ptr<Eigen::VectorXd> vectorPart, std::string label) {
  /*
   * Check cache capacity.
   *    If cache is full. Store on disk.
   *    Else store in _vectorPool
   */
  assert(vectorPart && "Pointer to vector segment is invalid!");
  if (_vectorPool.find(label) != _vectorPool.end()) {
    // entry already existing. Delete old entry.
    deleteVectorSegmentFromController(label);
  }
  _segmentSizes.insert(std::pair<std::string, int>(label, vectorPart->size()));
  if ((_nDoubles + vectorPart->size()) * 8 > _cacheLimit) {
    // Cache limit reached. Write copy on disk.
    saveVectorSegment(*vectorPart, label);
    // remove segment from cache.
    _vectorPool.insert(std::pair<std::string, std::shared_ptr<Eigen::VectorXd>>(label, nullptr));
  }
  else {
    // Cache limit is not reached yet! Save copy in cache.
    _vectorPool.insert(
        std::pair<std::string, std::shared_ptr<Eigen::VectorXd>>(label, std::make_shared<Eigen::VectorXd>(*vectorPart)));
    _nDoubles += vectorPart->size();
  }
}

std::shared_ptr<Eigen::VectorXd> VectorOnDiskStorageController::getVectorSegment(std::string label) {
  /*
   * Search cache.
   * If not already loaded, load it.
   */
  if (_vectorPool.find(label) == _vectorPool.end()) {
    throw std::runtime_error("There is no vector stored with label: " + label);
  }
  auto vectorSegment = _vectorPool.at(label);
  if (!vectorSegment) {
    // only nullptr stored. Load from disk.
    vectorSegment = loadVectorSegment(label);
  }
  assert(vectorSegment && "Loading the vector segment was not successful");
  return vectorSegment;
}

double VectorOnDiskStorageController::operator*(VectorOnDiskStorageController& rhs) {
  /*
   * Loop over all vector segments in both vectors.
   * Build scalar product of segments and add them up.
   */
  assert(this->size() == rhs.size() && "The vectors are not of the same size.");
  double scalarProduct = 0.0;
  for (const auto& it : _vectorPool) {
    const auto& label = it.first;
    const auto lhsSegment = *this->getVectorSegment(label);
    const auto rhsSegment = *rhs.getVectorSegment(label);
    assert(lhsSegment.size() == rhsSegment.size() && "Vector segments are not of the same sizes");
    scalarProduct += lhsSegment.cwiseProduct(rhsSegment).sum();
  }
  return scalarProduct;
}

std::shared_ptr<Eigen::VectorXd> VectorOnDiskStorageController::loadVectorSegment(std::string label) {
  assert(_vectorPool[label] == nullptr && "Vector already loaded!");
  assert(_files[label] && "There was no HDF5 file created!");
  auto vectorSegment = std::make_shared<Eigen::VectorXd>(_segmentSizes[label]);
  // if the HDF5 file is not opened yet. Open it.
  HDF5::Filepath name("tmp/" + label + "_" + _fileName);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(*_files[label], label);
  HDF5::load(*_files[label], label, *vectorSegment);
  return vectorSegment;
}

void VectorOnDiskStorageController::saveVectorSegment(Eigen::VectorXd& vectorSegment, std::string label) {
  // create directory if it does not exist yet.
  makePath("tmp");

  if (_files.find(label) == _files.end()) {
    // No entry available yet.
    try {
      _files.insert(std::pair<std::string, std::shared_ptr<H5::H5File>>(
          label, std::make_shared<H5::H5File>(("tmp/" + label + "_" + _fileName).c_str(), H5F_ACC_TRUNC)));
    }
    catch (H5::FileIException&) {
      throw std::runtime_error("A HDF5 file for the system is already open or has not been closed correctly!");
    }
  }
  else {
    if (_files[label]) {
      // Delete old file.
      try {
        _files[label]->close();
      }
      catch (...) {
      }
      std::remove(("tmp/" + label + "_" + _fileName).c_str());
    }
    _files[label] = std::make_shared<H5::H5File>(("tmp/" + label + "_" + _fileName).c_str(), H5F_ACC_TRUNC);
  }
  HDF5::save(*_files[label], label, vectorSegment);
}

void VectorOnDiskStorageController::removeVectorSegmentFromCache(std::string label) {
  saveVectorSegment(*_vectorPool[label], label);
  _nDoubles -= _segmentSizes[label];
  _vectorPool[label] = nullptr;
}

void VectorOnDiskStorageController::deleteVectorSegmentFromController(std::string label) {
  if (_vectorPool[label])
    _nDoubles -= _segmentSizes[label];
  _vectorPool.erase(label);
  _segmentSizes.erase(label);
  // delete the file which holds the data.
  if (_files[label]) {
    try {
      _files[label]->close();
    }
    catch (...) {
    }
    std::remove(("tmp/" + label + "_" + _fileName).c_str());
  }
}

VectorOnDiskStorageController::~VectorOnDiskStorageController() {
  // Remove files and tmp directory.
  for (const auto& file : _files) {
    if (file.second) {
      try {
        file.second->close();
      }
      catch (...) {
      }
      std::remove(("tmp/" + file.first + "_" + _fileName).c_str());
    }
  }
  std::remove("tmp");
}

} /* namespace Serenity */
