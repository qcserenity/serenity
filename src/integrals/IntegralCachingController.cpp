/**
 * @file   IntegralCachingController.cpp
 *
 * @date   Sep 24, 2021
 * @author Moritz Bensberg, Johannes Scheffler
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
#include "integrals/IntegralCachingController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "basis/Shell.h"
#include "memory/MemoryManager.h"
/* Include Std and External Headers */
#if __linux__ || __unix__ || __unix
#include <malloc.h> //Free unused memory.
#endif
#include <vector>

namespace Serenity {

IntegralCachingController::IntegralCachingController(std::shared_ptr<BasisController> basisController, int intCondition)
  : _basisController(basisController),
    _intCondition(intCondition),
    _nThreads(omp_get_max_threads()),
    _memManager(MemoryManager::getInstance()),
    _memoryPerThread(_nThreads) {
  _basisController->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
}

// Cache an integral
void IntegralCachingController::cacheIntegral(unsigned ij, const double* intptr, size_t size, unsigned threadId) {
  (*_cacheVector)[ij].push_back(std::vector<double>{intptr, intptr + size});
  _memoryPerThread[threadId] -= sizeof(double) * size;
}

const double* IntegralCachingController::getIntegral(unsigned ij, size_t cdcounter) {
  if ((*_cacheVector)[ij].size() <= cdcounter) {
    return nullptr;
  }
  else {
    return std::move(&((*_cacheVector)[ij][cdcounter][0]));
  }
}

bool IntegralCachingController::checkMem() {
  if (_memAvailable) {
    double memorySum = 0;
    for (unsigned i = 0; i < _nThreads; i++) {
      memorySum += _memoryPerThread[i];
    }
    if (memorySum < 0) {
      _memAvailable = false;
    }
  }
  return _memAvailable;
}

bool IntegralCachingController::timeCondition(const Shell& aa, const Shell& bb, const Shell& cc, const Shell& dd) {
  int nA = aa.getNPrimitives();
  int nB = bb.getNPrimitives();
  int nC = cc.getNPrimitives();
  int nD = dd.getNPrimitives();

  int lA = aa.getAngularMomentum();
  int lB = bb.getAngularMomentum();
  int lC = cc.getAngularMomentum();
  int lD = dd.getAngularMomentum();

  int tl = 1;
  int sum = lA + lB + lC + lD;
  if (sum == 0) {
    tl = 6;
  }
  else if (sum == 1) {
    tl = 3;
  }
  else if (sum == 2) {
    tl = 2;
  }
  return (nA * nB * nC * nD * tl >= _intCondition);
}

void IntegralCachingController::creatMemManagingVec() {
  // Reserve 500 MB for each thread plus some extra space.
  double freeMem = (double)_memManager->getAvailableSystemMemory() * 0.8 - 5e+8 * _nThreads;
  double freeMemPerThread = freeMem / _nThreads;
  if (freeMem > 0) {
    for (unsigned iThread = 0; iThread < _nThreads; ++iThread) {
      _memoryPerThread[iThread] = freeMemPerThread;
    }
  }
  else {
    for (unsigned iThread = 0; iThread < _nThreads; ++iThread) {
      _memoryPerThread[iThread] = -1;
    }
  }
}

void IntegralCachingController::clearCache() {
  _cacheVector = nullptr;
  // Free unused memory.
#if __linux__ || __unix__ || __unix
  malloc_trim(0);
#endif
}

void IntegralCachingController::notify() {
  unsigned nShells = _basisController->getReducedNBasisFunctions();
  this->clearCache();
  _cacheVector = std::make_unique<std::vector<std::vector<std::vector<double>>>>(nShells * (nShells + 1) / 2);
}
void IntegralCachingController::setPrescreeningThreshold(double threshold) {
  if (_prescreeningThreshold != threshold || !_cacheVector) {
    notify();
  }
  _prescreeningThreshold = threshold;
}

IntegralCachingController::~IntegralCachingController() {
  clearCache();
}

} /* namespace Serenity */
