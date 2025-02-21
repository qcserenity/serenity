/**
 * @file TwoElecThreeCenterCalculator.cpp
 *
 * @date Jul 12, 2021
 * @author Niklas Niemeyer
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
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
/* Include Serenity Internal Headers */
#include "memory/MemoryManager.h"

namespace Serenity {
TwoElecThreeCenterCalculator::TwoElecThreeCenterCalculator(LIBINT_OPERATOR op, double mu,
                                                           std::shared_ptr<BasisController> basisA,
                                                           std::shared_ptr<BasisController> basisB,
                                                           std::shared_ptr<BasisController> auxbasis,
                                                           double prescreeningThreshold, double maxDens)
  : _op(op),
    _basisControllerA(basisA),
    _basisControllerB(basisB == nullptr ? basisA : basisB),
    _twoBasisMode(_basisControllerA != _basisControllerB),
    _auxbasis(auxbasis),
    _prescreeningThreshold(prescreeningThreshold),
    _mu(mu),
    _nb_A(_basisControllerA->getNBasisFunctions()),
    _nb_B(_basisControllerB->getNBasisFunctions()),
    _nss(0),
    _offsets(auxbasis->getReducedNBasisFunctions(), 0) {
  _integrals.resize(omp_get_max_threads());

  // Initialize Libint for this calculator to use.
  unsigned maxPrimA = _basisControllerA->getMaxNumberOfPrimitives();
  unsigned maxPrimB = _basisControllerB->getMaxNumberOfPrimitives();
  unsigned maxPrimC = _auxbasis->getMaxNumberOfPrimitives();
  unsigned maxPrim = std::max({maxPrimA, maxPrimB, maxPrimC});
  _libint->initialize(_op, 0, 3, std::vector<std::shared_ptr<Atom>>(0), mu, std::numeric_limits<double>::epsilon(),
                      maxDens, maxPrim);

  // Setup shell pair data.
  this->setupShellPairs();

  // This is important to ensure that no race conditions occur when calculateIntegrals()
  // is called by multiple threads and the prescreening factors are not yet present.
  _auxbasis->getRIPrescreeningFactors();
}

TwoElecThreeCenterCalculator::~TwoElecThreeCenterCalculator() {
  _libint->finalize(_op, 0, 3);
}

Eigen::Ref<Eigen::MatrixXd> TwoElecThreeCenterCalculator::calculateIntegrals(unsigned P, unsigned iThread) {
  // Get basis information.
  auto& basisA = _basisControllerA->getBasis();
  auto& basisB = _basisControllerB->getBasis();
  auto& auxbasis = _auxbasis->getBasis();
  bool normAux = !(_auxbasis->isAtomicCholesky());

  const auto& auxPrescreening = _auxbasis->getRIPrescreeningFactors();
  unsigned np = auxbasis[P]->getNContracted();

  // Prepare integral loop.
  Eigen::MatrixXd dummy;
  _integrals[iThread] = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(_nb_A * _nb_B, np));
  double* matptr = _integrals[iThread]->data();

  for (auto& MN : (*_shellPairData)) {
    if (MN.factor * (*auxPrescreening)[P].factor < _prescreeningThreshold) {
      break;
    }

    // Shell indices of regular basis.
    unsigned I = MN.bf1;
    unsigned J = MN.bf2;

    // Number of basis functions in I and J shells.
    unsigned ni = basisA[I]->getNContracted();
    unsigned nj = basisB[J]->getNContracted();
    if (_libint->compute(_op, 0, *auxbasis[P], *basisB[J], *basisA[I], dummy, normAux)) {
      const double* intptr = dummy.data();
      for (unsigned p = 0; p < np; ++p) {
        unsigned long offset = p * _nb_A * _nb_B;
        for (unsigned j = 0; j < nj; ++j) {
          unsigned jj = _basisControllerB->extendedIndex(J) + j;
          for (unsigned i = 0; i < ni; ++i, ++intptr) {
            unsigned ii = _basisControllerA->extendedIndex(I) + i;
            matptr[ii + jj * _nb_A + offset] = (*intptr);
            if (!_twoBasisMode) {
              matptr[jj + ii * _nb_A + offset] = (*intptr);
            }
          }
        }
      }
    }
  }

  return *_integrals[iThread];
}

void TwoElecThreeCenterCalculator::setupShellPairs() {
  if (_twoBasisMode) {
    double maxPrim = std::max(_basisControllerA->getMaxNumberOfPrimitives(), _basisControllerB->getMaxNumberOfPrimitives());
    _libint->initialize(_op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), _mu, std::numeric_limits<double>::epsilon(),
                        10, maxPrim);
    _shellPairData = std::make_shared<std::vector<ShellPairData>>();
    auto& basisA = _basisControllerA->getBasis();
    auto& basisB = _basisControllerB->getBasis();
    Eigen::MatrixXd integrals;
    for (unsigned I = 0; I < basisA.size(); ++I) {
      auto& shellI = *basisA[I];
      for (unsigned J = 0; J < basisB.size(); ++J) {
        auto& shellJ = *basisB[J];
        if (_libint->compute(_op, 0, shellI, shellJ, shellI, shellJ, integrals)) {
          _shellPairData->push_back(ShellPairData(I, J, std::sqrt(integrals.maxCoeff()), false));
        }
      }
    }
    std::sort(_shellPairData->begin(), _shellPairData->end());
    std::reverse(_shellPairData->begin(), _shellPairData->end());
    _libint->finalize(_op, 0, 4);
  }
  else {
    _shellPairData = _basisControllerA->getShellPairData();
  }
}

void TwoElecThreeCenterCalculator::cacheIntegrals() {
  auto mem = MemoryManager::getInstance();
  double freeMem = 0.5 * mem->getAvailableSystemMemory();
  auto& auxbasis = _auxbasis->getBasis();

  double memx = sizeof(double) * _nb_A * _nb_B;
  size_t nxs = std::round(std::min((double)_auxbasis->getNBasisFunctions(), freeMem / memx));
  printf("  Caching %5i aux. basis functions (%3.0f%%, %5.2f GB).\n\n", (int)nxs,
         100 * (double)nxs / (double)_auxbasis->getNBasisFunctions(), 1e-9 * memx * nxs);

  // Determine number of shells stored.
  size_t stop = 0;
  for (size_t iShell = 0; iShell < _auxbasis->getReducedNBasisFunctions(); ++iShell) {
    if (iShell + 1 < _offsets.size()) {
      _offsets[iShell + 1] = _offsets[iShell] + auxbasis[iShell]->getNContracted();
    }
    if (stop + auxbasis[iShell]->getNContracted() > nxs) {
      break;
    }
    stop += auxbasis[iShell]->getNContracted();
    _nss += 1;
  }

  _cache.resize(_nb_A * _nb_B, nxs);

#pragma omp parallel for schedule(dynamic)
  for (size_t iShell = 0; iShell < _nss; ++iShell) {
    auto integrals = this->calculateIntegrals(iShell, omp_get_thread_num());
    _cache.middleCols(_offsets[iShell], integrals.cols()) = integrals;
  }
}

void TwoElecThreeCenterCalculator::clearCache() {
  _nss = 0;
  _cache.resize(0, 0);
}

} /* namespace Serenity */
