/**
 * @file OrbitalPairSet.cpp
 *
 * @date Feb. 18, 2021
 * @author Moritz Bensberg
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
#include "data/OrbitalPairSet.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"   //Load/Write integrals.
#include "io/HDF5.h"            //Write to disk.
#include "misc/SerenityError.h" //Errors.
/* Include Std and External Headers */
#if __linux__ || __unix__ || __unix
#include <malloc.h> //Free unused memory.
#endif
namespace Serenity {

OrbitalPairSet::OrbitalPairSet() = default;

OrbitalPairSet::~OrbitalPairSet() = default;

void OrbitalPairSet::toHDF5(HDF5::H5File& file) {
  _onDisk = true;
  for (auto& pair : *this) {
    pair->writeIntegralsToFile(file);
  }
  file.flush(H5F_SCOPE_GLOBAL);
}

void OrbitalPairSet::fromHDF5(HDF5::H5File& file) {
  if (!_onDisk)
    throw SerenityError("ERROR: Data was never written to disk. It is impossible to read it!");
  for (auto& pair : *this) {
    pair->loadIntegralsFromFile(file);
  }
  _inMemory = true;
}

void OrbitalPairSet::removeInteralsFromMemory() {
  for (auto& pair : *this) {
    pair->flushIntegrals();
  }
  _inMemory = false;
  // Free unused memory.
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

void OrbitalPairSet::setInMemory(bool inMemory) {
  _inMemory = inMemory;
}

bool OrbitalPairSet::integralsReady() {
  return _inMemory;
}

double OrbitalPairSet::memoryDemand(bool sigmaInts) {
  if (!_memory) {
    _memory = std::make_shared<double>(0.0);
    for (auto& pair : *this)
      *_memory += pair->getMemoryRequirement(false, sigmaInts);
  }
  return *_memory;
}

Eigen::SparseVector<int> OrbitalPairSet::getTotalFittingDomain() {
  Eigen::SparseVector<int> totalKDomain(this->at(0)->getFittingDomain().rows());
  for (auto& pair : *this) {
    totalKDomain += pair->getFittingDomain();
  }
  return totalKDomain;
}

} /* namespace Serenity */