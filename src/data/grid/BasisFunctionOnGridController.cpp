/**
 * @file   BasisFunctionOnGridController.cpp
 *
 * @date   May 7, 2014, last rework June 6, 2017
 * @author Thomas Dresselhaus, last rework Jan Unsleber
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
#include "data/grid/BasisFunctionOnGridController.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "geometry/Point.h"
#include "grid/GridController.h"
#include "math/Derivatives.h"
#include "math/FloatMaths.h"
#include "misc/HelperFunctions.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <omp.h>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Serenity {

BasisFunctionOnGridController::BasisFunctionBlockOnGridData::BasisFunctionBlockOnGridData(const unsigned int nBasisFunctions,
                                                                                          const unsigned int blockSize,
                                                                                          const unsigned int derivativeLevel)
  : negligible(nBasisFunctions),
    functionValues(blockSize, nBasisFunctions),
    derivativeValues(derivativeLevel >= 1 ? makeGradientPtr<Eigen::MatrixXd>(blockSize, nBasisFunctions) : nullptr),
    secondDerivativeValues(derivativeLevel >= 2 ? makeHessianPtr<Eigen::MatrixXd>(blockSize, nBasisFunctions) : nullptr) {
  if (derivativeLevel > 2)
    throw SerenityError("Third derivatives of basis functions on the grid are not implemented!");
  negligible.setZero();
  functionValues.setZero();
  if (derivativeValues != nullptr) {
    derivativeValues->x.setZero();
    derivativeValues->y.setZero();
    derivativeValues->z.setZero();
  }
  if (secondDerivativeValues != nullptr) {
    secondDerivativeValues->xx.setZero();
    secondDerivativeValues->xy.setZero();
    secondDerivativeValues->xz.setZero();
    secondDerivativeValues->yy.setZero();
    secondDerivativeValues->yz.setZero();
    secondDerivativeValues->zz.setZero();
  }
}

BasisFunctionOnGridController::BasisFunctionBlockOnGridData::BasisFunctionBlockOnGridData(const BasisFunctionBlockOnGridData& orig)
  : center(orig.center),
    negligible(orig.negligible),
    functionValues(orig.functionValues),
    derivativeValues(orig.derivativeValues ? new Gradient<Eigen::MatrixXd>(*orig.derivativeValues) : nullptr),
    secondDerivativeValues(orig.secondDerivativeValues ? new Hessian<Eigen::MatrixXd>(*orig.secondDerivativeValues) : nullptr) {
}

BasisFunctionOnGridController::BasisFunctionOnGridController(std::shared_ptr<BasisController> basisController,
                                                             std::shared_ptr<GridController> gridController,
                                                             const unsigned int maxBlockSize, const double radialThreshold,
                                                             const unsigned int highestDerivative)
  : _basisController(basisController),
    _gridController(gridController),
    _maxBlockSize(maxBlockSize),
    _nPoints(_gridController->getNGridPoints()),
    _nBlocks((unsigned int)ceil((double)_nPoints / _maxBlockSize)),
    _radialThreshold(radialThreshold),
    _exponentThreshold(-log(_radialThreshold)),
    _highestDerivative(highestDerivative),
    _lastBlock(nullptr) {
  if (highestDerivative > 2)
    throw SerenityError("Third derivatives of basis functions on the grid are not implemented!");
#ifdef _OPENMP
  _workspace.resize(omp_get_max_threads());
  for (int i = 0; i < omp_get_max_threads(); ++i) {
#else
  _workspace.resize(1);
  for (unsigned int i = 0; i < 1; ++i) {
#endif
    _workspace[i].reset(
        new BasisFunctionBlockOnGridData(_basisController->getNBasisFunctions(), _maxBlockSize, _highestDerivative));
  }
  _gridController->addSensitiveObject(this->_self);
}

void BasisFunctionOnGridController::notify() {
  _upToDate = false;
}

void BasisFunctionOnGridController::setHighestDerivative(unsigned int newHighestDerivative) {
  if (newHighestDerivative > 2)
    throw SerenityError("Third derivatives of basis functions on the grid are not implemented!");
  if (newHighestDerivative > _highestDerivative) {
    std::lock_guard<std::mutex> lock(_lock);
    /* Reinit workspace */
    _workspace.resize(omp_get_max_threads());
    for (int i = 0; i < omp_get_max_threads(); ++i) {
      _workspace[i].reset(
          new BasisFunctionBlockOnGridData(_basisController->getNBasisFunctions(), _maxBlockSize, newHighestDerivative));
    }
  }
  _highestDerivative = newHighestDerivative;
}

const std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData>
BasisFunctionOnGridController::getBlockOnGridData(const unsigned int blockIndex) {
  if (!_upToDate) {
    _nPoints = _gridController->getNGridPoints();
    _nBlocks = (unsigned int)ceil((double)_nPoints / _maxBlockSize);
    _upToDate = true;
  }
  assert(blockIndex < _nBlocks);
  return calculateBasisFunctionData(blockIndex);
}

unsigned int BasisFunctionOnGridController::getFirstIndexOfBlock(const unsigned int blockIndex) {
  if (!_upToDate) {
    _nPoints = _gridController->getNGridPoints();
    _nBlocks = (unsigned int)ceil((double)_nPoints / _maxBlockSize);
    _upToDate = true;
  }
  return blockIndex * _maxBlockSize;
}

unsigned int BasisFunctionOnGridController::getNGridPoints() const {
  return _gridController->getNGridPoints();
}

unsigned int BasisFunctionOnGridController::getNBasisFunctions() const {
  return _basisController->getNBasisFunctions();
}

std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData>
BasisFunctionOnGridController::calculateBasisFunctionData(const unsigned int blockIndex) {
  if (!_upToDate) {
    _nPoints = _gridController->getNGridPoints();
    _nBlocks = (unsigned int)ceil((double)_nPoints / _maxBlockSize);
    _upToDate = true;
  }
  assert(blockIndex < _nBlocks);

  if (omp_get_thread_num() == 0)
    Timings::takeTime("Tech. -    Basis On Grid Eval.");

  /*
   * Get in data locally
   */
  const unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
  const unsigned int nReducedBasisFunctions = _basisController->getReducedNBasisFunctions();
  auto& gridPoints = _gridController->getGridPoints();
  const double* grid_ptr = gridPoints.data();
  const unsigned int firstIndex = getFirstIndexOfBlock(blockIndex);
  unsigned int blockSize;
  if (blockIndex == _nBlocks - 1) {
    blockSize = _nPoints % _maxBlockSize;
    if (blockSize == 0) {
      blockSize = _maxBlockSize;
    }
  }
  else {
    blockSize = _maxBlockSize;
  }
  const auto& basis = _basisController->getBasis();
  /*
   * Create object to store data
   */
  std::shared_ptr<BasisFunctionBlockOnGridData> thisBlockData;
  if (blockIndex == _nBlocks - 1) {
    _lastBlock.reset(new BasisFunctionBlockOnGridData(nBasisFunctions, blockSize, _highestDerivative));
    thisBlockData = _lastBlock;
  }
  else {
    thisBlockData = _workspace[omp_get_thread_num()];
    assert(blockSize == _maxBlockSize);
  }
  double* values = thisBlockData->functionValues.data();
  Gradient<double*> derivatives;
  if (_highestDerivative >= 1)
    zip(*thisBlockData->derivativeValues, derivatives, [](Eigen::MatrixXd& matrix, double*& ptr) { ptr = matrix.data(); });
  Hessian<double*> secondDerivatives;
  if (_highestDerivative >= 2)
    zip(*thisBlockData->secondDerivativeValues, secondDerivatives,
        [](Eigen::MatrixXd& matrix, double*& ptr) { ptr = matrix.data(); });

  /*
   * Loop over basis functions and prescreen
   *
   * Here the entire block is prescreened
   * If a basis function is not significant on the center
   * of this block, it will not be calculated at all for this block.
   * It is deemed 'negligible'.
   *
   */
  Eigen::VectorXd center = _gridController->getGridPoints().block(0, firstIndex, 3, blockSize).rowwise().mean();
  const double* c_ptr = center.data();
  thisBlockData->center = center;
  double spread =
      (_gridController->getGridPoints().block(0, firstIndex, 3, blockSize).colwise() - center).colwise().norm().maxCoeff();
  thisBlockData->negligible.setZero();
  for (unsigned int muRed = 0; muRed < nReducedBasisFunctions; muRed++) {
    // Gather data
    const auto& basMu = basis[muRed];
    const unsigned int nPrimMu = basMu->getNPrimitives();
    const double* expMu = basMu->alpha.data();
    const auto& primNormsMu = basMu->contr[0].coeff;

    // calculate distance between bf center and block center
    const double pMuX(c_ptr[0] - basMu->getX());
    const double pMuY(c_ptr[1] - basMu->getY());
    const double pMuZ(c_ptr[2] - basMu->getZ());
    double distPMuSquare = sqrt(pMuX * pMuX + pMuY * pMuY + pMuZ * pMuZ) - spread;

    /*
     * If the sphere around the center of the gridblock created by the "spread" is within
     * 1 au of the basis function (BF) center, it is assumed that the BF is not negligible.
     * Thus, the rest of the negligibility check is skipped.
     */
    if (distPMuSquare < 1.0)
      continue;
    distPMuSquare = distPMuSquare * distPMuSquare;

    // evaluate radial part of this shell
    double radial = 0.0;
    for (unsigned int i = 0; i < nPrimMu; i++) {
      const double tmp = expMu[i] * distPMuSquare;
      const double expMinusTmp = exp(-tmp);
      radial += primNormsMu[i] * expMinusTmp;
    }

    // set all bfs in shell negligible if needed
    if (fabs(radial) < _radialThreshold) {
      unsigned int start = _basisController->extendedIndex(muRed);
      unsigned int stop = start + basMu->getNContracted();
      for (unsigned int i = start; i < stop; i++) {
        thisBlockData->negligible[i] = 1;
      }
    }
  }
  /*
   * Loop over points in this block
   */
  for (unsigned int iBlock = 0; iBlock < blockSize; iBlock++) {
    const double& pX = grid_ptr[(firstIndex + iBlock) * 3 + 0];
    const double& pY = grid_ptr[(firstIndex + iBlock) * 3 + 1];
    const double& pZ = grid_ptr[(firstIndex + iBlock) * 3 + 2];

    /*
     * Loop over basis functions and calculate if not negligible
     */
    for (unsigned int muRed = 0; muRed < nReducedBasisFunctions; muRed++) {
      if (thisBlockData->negligible[_basisController->extendedIndex(muRed)])
        continue;

      /*
       * Get in data locally
       */
      const auto& basMu = basis[muRed];
      // number of different angular parts for Phi_mu
      const unsigned int nMu = basMu->getNContracted();
      const unsigned int nPrimMu = basMu->getNPrimitives();
      const auto& primNormsMu = basMu->contr[0].coeff;
      const auto& normFactorOfMu = basMu->getNormFactors();
      const double* expMu = basMu->alpha.data();

      const unsigned int angularMomentumOfMu = basMu->getAngularMomentum();
      const unsigned int iMuStart = _basisController->extendedIndex(muRed);
      /*
       * Calc Phi_mu(p)
       */
      /*
       * Radial part
       */
      const double pMuX(pX - basMu->getX());
      const double pMuY(pY - basMu->getY());
      const double pMuZ(pZ - basMu->getZ());
      const double distPMuSquare = pMuX * pMuX + pMuY * pMuY + pMuZ * pMuZ;

      double radial = 0.0;
      double dradial = 0.0;
      double ddradial = 0.0;
      for (unsigned int i = 0; i < nPrimMu; i++) {
        const double tmp = expMu[i] * distPMuSquare;
        if (tmp < _exponentThreshold) {
          const double expMinusTmp = exp(-tmp);
          radial += primNormsMu[i] * expMinusTmp;
          dradial -= 2.0 * expMu[i] * primNormsMu[i] * expMinusTmp;
          ddradial += 4.0 * expMu[i] * expMu[i] * primNormsMu[i] * expMinusTmp;
        }
      }
      // The first index in phiOfMuAndP for this shell and grid point
      const unsigned int iMuPStart = iMuStart * blockSize + iBlock;
      /*
       * Prescreen function on this point
       */
      if (fabs(radial) < _radialThreshold) {
        /*
         * This shell / point combination is insignificant. Skip the angular part.
         */
        const unsigned int iMuPEnd = iMuPStart + nMu * blockSize;
        for (unsigned int iMuP = iMuPStart; iMuP < iMuPEnd; iMuP += blockSize) {
          values[iMuP] = 0.0;
          if (_highestDerivative >= 1) {
            for (auto& component : derivatives)
              component[iMuP] = 0.0;
            if (_highestDerivative >= 2) {
              for (auto& component : secondDerivatives)
                component[iMuP] = 0.0;
            }
          }
        }
        continue;
      }
      /*
       * Angular part
       */
      double x[AM_MAX + 3];
      double y[AM_MAX + 3];
      double z[AM_MAX + 3];
      x[0] = 1.0;
      y[0] = 1.0;
      z[0] = 1.0;
      for (unsigned int exponent = 1; exponent <= angularMomentumOfMu + _highestDerivative; exponent++) {
        x[exponent] = x[exponent - 1] * pMuX;
        y[exponent] = y[exponent - 1] * pMuY;
        z[exponent] = z[exponent - 1] * pMuZ;
      }
      double* valueStore = nullptr;
      auto derivativeStore = makeGradient<double*>(nullptr);
      auto hessianStore = makeHessian<double*>(nullptr);
      unsigned int offset;

      // Store the values directly
      valueStore = values + iMuPStart;
      if (_highestDerivative >= 1) {
        zip(derivativeStore, derivatives, [&iMuPStart](double*& store, double*& orig) { store = orig + iMuPStart; });
        if (_highestDerivative >= 2) {
          zip(hessianStore, secondDerivatives, [&iMuPStart](double*& store, double*& orig) { store = orig + iMuPStart; });
        }
      }
      offset = blockSize;

      if (!basMu->isSpherical()) {
        for (int muX = angularMomentumOfMu, iMuP = 0, iMu = iMuStart, iMuShell = 0; muX >= 0; muX--) {
          for (int muY = angularMomentumOfMu - muX; muY >= 0; muY--, iMuP += offset, iMu++, iMuShell++) {
            const int muZ = angularMomentumOfMu - muY - muX;

            const double value = x[muX] * y[muY] * z[muZ] * radial * normFactorOfMu[iMuShell];

            if (_highestDerivative >= 1) {
              derivativeStore.x[iMuP] = dradial * x[muX + 1] * y[muY] * z[muZ] * normFactorOfMu[iMuShell];
              if (muX > 0) {
                derivativeStore.x[iMuP] += muX * x[muX - 1] * y[muY] * z[muZ] * radial * normFactorOfMu[iMuShell];
              }
              derivativeStore.y[iMuP] = dradial * x[muX] * y[muY + 1] * z[muZ] * normFactorOfMu[iMuShell];
              if (muY > 0) {
                derivativeStore.y[iMuP] += muY * x[muX] * y[muY - 1] * z[muZ] * radial * normFactorOfMu[iMuShell];
              }
              derivativeStore.z[iMuP] = dradial * x[muX] * y[muY] * z[muZ + 1] * normFactorOfMu[iMuShell];
              if (muZ > 0) {
                derivativeStore.z[iMuP] += muZ * x[muX] * y[muY] * z[muZ - 1] * radial * normFactorOfMu[iMuShell];
              }

              if (_highestDerivative >= 2) {
                // XX
                hessianStore.xx[iMuP] = (ddradial * x[muX + 2] + dradial * x[muX] * (2 * muX + 1)) * y[muY] * z[muZ] *
                                        normFactorOfMu[iMuShell];
                if (muX > 1) {
                  hessianStore.xx[iMuP] += muX * (muX - 1) * x[muX - 2] * y[muY] * z[muZ] * radial * normFactorOfMu[iMuShell];
                }
                // XY
                hessianStore.xy[iMuP] = ddradial * x[muX + 1] * y[muY + 1] * z[muZ] * normFactorOfMu[iMuShell];
                if (muX > 0) {
                  hessianStore.xy[iMuP] += dradial * muX * x[muX - 1] * y[muY + 1] * z[muZ] * normFactorOfMu[iMuShell];
                }
                if (muY > 0) {
                  hessianStore.xy[iMuP] += dradial * muY * x[muX + 1] * y[muY - 1] * z[muZ] * normFactorOfMu[iMuShell];
                  if (muX > 0) {
                    hessianStore.xy[iMuP] += muX * muY * x[muX - 1] * y[muY - 1] * z[muZ] * normFactorOfMu[iMuShell] * radial;
                  }
                }
                // XZ
                hessianStore.xz[iMuP] = x[muX + 1] * y[muY] * z[muZ + 1] * ddradial * normFactorOfMu[iMuShell];
                if (muX > 0) {
                  hessianStore.xz[iMuP] += dradial * muX * x[muX - 1] * y[muY] * z[muZ + 1] * normFactorOfMu[iMuShell];
                }
                if (muZ > 0) {
                  hessianStore.xz[iMuP] += dradial * muZ * x[muX + 1] * y[muY] * z[muZ - 1] * normFactorOfMu[iMuShell];
                  if (muX > 0) {
                    hessianStore.xz[iMuP] += muX * muZ * x[muX - 1] * y[muY] * z[muZ - 1] * normFactorOfMu[iMuShell] * radial;
                  }
                }
                // YY
                hessianStore.yy[iMuP] = (ddradial * y[muY + 2] + dradial * y[muY] * (2 * muY + 1)) * x[muX] * z[muZ] *
                                        normFactorOfMu[iMuShell];
                if (muY > 1) {
                  hessianStore.yy[iMuP] +=

                      muY * (muY - 1) * y[muY - 2] * x[muX] * z[muZ] * radial * normFactorOfMu[iMuShell];
                }
                // YZ
                hessianStore.yz[iMuP] = ddradial * x[muX] * y[muY + 1] * z[muZ + 1] * normFactorOfMu[iMuShell];
                if (muY > 0) {
                  hessianStore.yz[iMuP] += dradial * muY * x[muX] * y[muY - 1] * z[muZ + 1] * normFactorOfMu[iMuShell];
                }
                if (muZ > 0) {
                  hessianStore.yz[iMuP] += dradial * muZ * x[muX] * y[muY + 1] * z[muZ - 1] * normFactorOfMu[iMuShell];
                  if (muY > 0) {
                    hessianStore.yz[iMuP] += muY * muZ * x[muX] * y[muY - 1] * z[muZ - 1] * normFactorOfMu[iMuShell] * radial;
                  }
                }
                // ZZ
                hessianStore.zz[iMuP] = (ddradial * z[muZ + 2] + dradial * z[muZ] * (2 * muZ + 1)) * y[muY] * x[muX] *
                                        normFactorOfMu[iMuShell];
                if (muZ > 1) {
                  hessianStore.zz[iMuP] += muZ * (muZ - 1) * z[muZ - 2] * y[muY] * x[muX] * radial * normFactorOfMu[iMuShell];
                }
              }

            } /**/
            valueStore[iMuP] = value;
          }
        }
      }
      else {
        /*
         * The following code for spherical basis functions
         * has been generated using Python3 code.
         * The code should reside in the same folder as this class
         * and is called sharmonics.py
         * The code uses the Herglotz generating function including separation
         * into z+r and x+y factors for real harmonics in Cartesian space.
         * - JU
         */

        /* ==================
         *   Initialization
         * ==================*/
        double Y[2 * AM_MAX + 1];
        double dYdx[2 * AM_MAX + 1];
        double dYdy[2 * AM_MAX + 1];
        double dYdz[2 * AM_MAX + 1];
        double d2Ydxdx[2 * AM_MAX + 1];
        double d2Ydxdy[2 * AM_MAX + 1];
        double d2Ydxdz[2 * AM_MAX + 1];
        double d2Ydydy[2 * AM_MAX + 1];
        double d2Ydydz[2 * AM_MAX + 1];
        double d2Ydzdz[2 * AM_MAX + 1];

        /* =============
         *   Harmonics
         * =============*/
        switch (angularMomentumOfMu) {
          case 0:
            Y[0] = 1.0;
            if (_highestDerivative >= 1) {
              dYdx[0] = 0.0;
              dYdy[0] = 0.0;
              dYdz[0] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = 0.0;
                d2Ydxdy[0] = 0.0;
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = 0.0;
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
              }
            }
            break;
          case 1:
            Y[0] = y[1];
            Y[1] = z[1];
            Y[2] = x[1];
            if (_highestDerivative >= 1) {
              dYdx[0] = 0.0;
              dYdy[0] = 1.0;
              dYdz[0] = 0.0;
              dYdx[1] = 0.0;
              dYdy[1] = 0.0;
              dYdz[1] = 1.0;
              dYdx[2] = 1.0;
              dYdy[2] = 0.0;
              dYdz[2] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = 0.0;
                d2Ydxdy[0] = 0.0;
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = 0.0;
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = 0.0;
                d2Ydxdy[1] = 0.0;
                d2Ydxdz[1] = 0.0;
                d2Ydydy[1] = 0.0;
                d2Ydydz[1] = 0.0;
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = 0.0;
                d2Ydxdy[2] = 0.0;
                d2Ydxdz[2] = 0.0;
                d2Ydydy[2] = 0.0;
                d2Ydydz[2] = 0.0;
                d2Ydzdz[2] = 0.0;
              }
            }
            break;
          case 2:
            Y[0] = sqrt(3.0 / 4.0) * (2.0 * x[1] * y[1]);
            Y[1] = sqrt(3.0) * (1.0 * y[1] * z[1]);
            Y[2] = sqrt(1.0 / 4.0) * (2.0 * z[2] - 1.0 * x[2] - 1.0 * y[2]);
            Y[3] = sqrt(3.0) * (1.0 * x[1] * z[1]);
            Y[4] = sqrt(3.0 / 4.0) * (1.0 * x[2] - 1.0 * y[2]);
            if (_highestDerivative >= 1) {
              dYdx[0] = sqrt(3.0 / 4.0) * (2.0 * y[1]);
              dYdy[0] = sqrt(3.0 / 4.0) * (2.0 * x[1]);
              dYdz[0] = 0.0;
              dYdx[1] = 0.0;
              dYdy[1] = sqrt(3.0) * (1.0 * z[1]);
              dYdz[1] = sqrt(3.0) * (1.0 * y[1]);
              dYdx[2] = sqrt(1.0 / 4.0) * (-2.0 * x[1]);
              dYdy[2] = sqrt(1.0 / 4.0) * (-2.0 * y[1]);
              dYdz[2] = sqrt(1.0 / 4.0) * (4.0 * z[1]);
              dYdx[3] = sqrt(3.0) * (1.0 * z[1]);
              dYdy[3] = 0.0;
              dYdz[3] = sqrt(3.0) * (1.0 * x[1]);
              dYdx[4] = sqrt(3.0 / 4.0) * (2.0 * x[1]);
              dYdy[4] = sqrt(3.0 / 4.0) * (-2.0 * y[1]);
              dYdz[4] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = 0.0;
                d2Ydxdy[0] = sqrt(3.0 / 4.0) * (2.0);
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = 0.0;
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = 0.0;
                d2Ydxdy[1] = 0.0;
                d2Ydxdz[1] = 0.0;
                d2Ydydy[1] = 0.0;
                d2Ydydz[1] = sqrt(3.0) * (1.0);
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = sqrt(1.0 / 4.0) * (-2.0);
                d2Ydxdy[2] = 0.0;
                d2Ydxdz[2] = 0.0;
                d2Ydydy[2] = sqrt(1.0 / 4.0) * (-2.0);
                d2Ydydz[2] = 0.0;
                d2Ydzdz[2] = sqrt(1.0 / 4.0) * (4.0);
                d2Ydxdx[3] = 0.0;
                d2Ydxdy[3] = 0.0;
                d2Ydxdz[3] = sqrt(3.0) * (1.0);
                d2Ydydy[3] = 0.0;
                d2Ydydz[3] = 0.0;
                d2Ydzdz[3] = 0.0;
                d2Ydxdx[4] = sqrt(3.0 / 4.0) * (2.0);
                d2Ydxdy[4] = 0.0;
                d2Ydxdz[4] = 0.0;
                d2Ydydy[4] = sqrt(3.0 / 4.0) * (-2.0);
                d2Ydydz[4] = 0.0;
                d2Ydzdz[4] = 0.0;
              }
            }
            break;
          case 3:
            Y[0] = sqrt(5.0 / 8.0) * (3.0 * x[2] * y[1] - 1.0 * y[3]);
            Y[1] = sqrt(15.0 / 4.0) * (2.0 * x[1] * y[1] * z[1]);
            Y[2] = sqrt(3.0 / 8.0) * (4.0 * y[1] * z[2] - 1.0 * x[2] * y[1] - 1.0 * y[3]);
            Y[3] = sqrt(1.0 / 4.0) * (2.0 * z[3] - 3.0 * x[2] * z[1] - 3.0 * y[2] * z[1]);
            Y[4] = sqrt(3.0 / 8.0) * (4.0 * x[1] * z[2] - 1.0 * x[3] - 1.0 * x[1] * y[2]);
            Y[5] = sqrt(15.0 / 4.0) * (1.0 * x[2] * z[1] - 1.0 * y[2] * z[1]);
            Y[6] = sqrt(5.0 / 8.0) * (1.0 * x[3] - 3.0 * x[1] * y[2]);
            if (_highestDerivative >= 1) {
              dYdx[0] = sqrt(5.0 / 8.0) * (6.0 * x[1] * y[1]);
              dYdy[0] = sqrt(5.0 / 8.0) * (3.0 * x[2] - 3.0 * y[2]);
              dYdz[0] = 0.0;
              dYdx[1] = sqrt(15.0 / 4.0) * (2.0 * y[1] * z[1]);
              dYdy[1] = sqrt(15.0 / 4.0) * (2.0 * x[1] * z[1]);
              dYdz[1] = sqrt(15.0 / 4.0) * (2.0 * x[1] * y[1]);
              dYdx[2] = sqrt(3.0 / 8.0) * (-2.0 * x[1] * y[1]);
              dYdy[2] = sqrt(3.0 / 8.0) * (4.0 * z[2] - 1.0 * x[2] - 3.0 * y[2]);
              dYdz[2] = sqrt(3.0 / 8.0) * (8.0 * y[1] * z[1]);
              dYdx[3] = sqrt(1.0 / 4.0) * (-6.0 * x[1] * z[1]);
              dYdy[3] = sqrt(1.0 / 4.0) * (-6.0 * y[1] * z[1]);
              dYdz[3] = sqrt(1.0 / 4.0) * (6.0 * z[2] - 3.0 * x[2] - 3.0 * y[2]);
              dYdx[4] = sqrt(3.0 / 8.0) * (4.0 * z[2] - 3.0 * x[2] - 1.0 * y[2]);
              dYdy[4] = sqrt(3.0 / 8.0) * (-2.0 * x[1] * y[1]);
              dYdz[4] = sqrt(3.0 / 8.0) * (8.0 * x[1] * z[1]);
              dYdx[5] = sqrt(15.0 / 4.0) * (2.0 * x[1] * z[1]);
              dYdy[5] = sqrt(15.0 / 4.0) * (-2.0 * y[1] * z[1]);
              dYdz[5] = sqrt(15.0 / 4.0) * (1.0 * x[2] - 1.0 * y[2]);
              dYdx[6] = sqrt(5.0 / 8.0) * (3.0 * x[2] - 3.0 * y[2]);
              dYdy[6] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * y[1]);
              dYdz[6] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = sqrt(5.0 / 8.0) * (6.0 * y[1]);
                d2Ydxdy[0] = sqrt(5.0 / 8.0) * (6.0 * x[1]);
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = sqrt(5.0 / 8.0) * (-6.0 * y[1]);
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = 0.0;
                d2Ydxdy[1] = sqrt(15.0 / 4.0) * (2.0 * z[1]);
                d2Ydxdz[1] = sqrt(15.0 / 4.0) * (2.0 * y[1]);
                d2Ydydy[1] = 0.0;
                d2Ydydz[1] = sqrt(15.0 / 4.0) * (2.0 * x[1]);
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = sqrt(3.0 / 8.0) * (-2.0 * y[1]);
                d2Ydxdy[2] = sqrt(3.0 / 8.0) * (-2.0 * x[1]);
                d2Ydxdz[2] = 0.0;
                d2Ydydy[2] = sqrt(3.0 / 8.0) * (-6.0 * y[1]);
                d2Ydydz[2] = sqrt(3.0 / 8.0) * (8.0 * z[1]);
                d2Ydzdz[2] = sqrt(3.0 / 8.0) * (8.0 * y[1]);
                d2Ydxdx[3] = sqrt(1.0 / 4.0) * (-6.0 * z[1]);
                d2Ydxdy[3] = 0.0;
                d2Ydxdz[3] = sqrt(1.0 / 4.0) * (-6.0 * x[1]);
                d2Ydydy[3] = sqrt(1.0 / 4.0) * (-6.0 * z[1]);
                d2Ydydz[3] = sqrt(1.0 / 4.0) * (-6.0 * y[1]);
                d2Ydzdz[3] = sqrt(1.0 / 4.0) * (12.0 * z[1]);
                d2Ydxdx[4] = sqrt(3.0 / 8.0) * (-6.0 * x[1]);
                d2Ydxdy[4] = sqrt(3.0 / 8.0) * (-2.0 * y[1]);
                d2Ydxdz[4] = sqrt(3.0 / 8.0) * (8.0 * z[1]);
                d2Ydydy[4] = sqrt(3.0 / 8.0) * (-2.0 * x[1]);
                d2Ydydz[4] = 0.0;
                d2Ydzdz[4] = sqrt(3.0 / 8.0) * (8.0 * x[1]);
                d2Ydxdx[5] = sqrt(15.0 / 4.0) * (2.0 * z[1]);
                d2Ydxdy[5] = 0.0;
                d2Ydxdz[5] = sqrt(15.0 / 4.0) * (2.0 * x[1]);
                d2Ydydy[5] = sqrt(15.0 / 4.0) * (-2.0 * z[1]);
                d2Ydydz[5] = sqrt(15.0 / 4.0) * (-2.0 * y[1]);
                d2Ydzdz[5] = 0.0;
                d2Ydxdx[6] = sqrt(5.0 / 8.0) * (6.0 * x[1]);
                d2Ydxdy[6] = sqrt(5.0 / 8.0) * (-6.0 * y[1]);
                d2Ydxdz[6] = 0.0;
                d2Ydydy[6] = sqrt(5.0 / 8.0) * (-6.0 * x[1]);
                d2Ydydz[6] = 0.0;
                d2Ydzdz[6] = 0.0;
              }
            }
            break;
          case 4:
            Y[0] = sqrt(35.0 / 64.0) * (4.0 * x[3] * y[1] - 4.0 * x[1] * y[3]);
            Y[1] = sqrt(35.0 / 8.0) * (3.0 * x[2] * y[1] * z[1] - 1.0 * y[3] * z[1]);
            Y[2] = sqrt(5.0 / 16.0) * (12.0 * x[1] * y[1] * z[2] - 2.0 * x[3] * y[1] - 2.0 * x[1] * y[3]);
            Y[3] = sqrt(5.0 / 8.0) * (4.0 * y[1] * z[3] - 3.0 * x[2] * y[1] * z[1] - 3.0 * y[3] * z[1]);
            Y[4] = sqrt(1.0 / 64.0) *
                   (8.0 * z[4] - 24.0 * x[2] * z[2] - 24.0 * y[2] * z[2] + 3.0 * x[4] + 3.0 * y[4] + 6.0 * x[2] * y[2]);
            Y[5] = sqrt(5.0 / 8.0) * (4.0 * x[1] * z[3] - 3.0 * x[3] * z[1] - 3.0 * x[1] * y[2] * z[1]);
            Y[6] = sqrt(5.0 / 16.0) * (6.0 * x[2] * z[2] - 1.0 * x[4] - 6.0 * y[2] * z[2] + 1.0 * y[4]);
            Y[7] = sqrt(35.0 / 8.0) * (1.0 * x[3] * z[1] - 3.0 * x[1] * y[2] * z[1]);
            Y[8] = sqrt(35.0 / 64.0) * (1.0 * x[4] - 6.0 * x[2] * y[2] + 1.0 * y[4]);
            if (_highestDerivative >= 1) {
              dYdx[0] = sqrt(35.0 / 64.0) * (12.0 * x[2] * y[1] - 4.0 * y[3]);
              dYdy[0] = sqrt(35.0 / 64.0) * (4.0 * x[3] - 12.0 * x[1] * y[2]);
              dYdz[0] = 0.0;
              dYdx[1] = sqrt(35.0 / 8.0) * (6.0 * x[1] * y[1] * z[1]);
              dYdy[1] = sqrt(35.0 / 8.0) * (3.0 * x[2] * z[1] - 3.0 * y[2] * z[1]);
              dYdz[1] = sqrt(35.0 / 8.0) * (3.0 * x[2] * y[1] - 1.0 * y[3]);
              dYdx[2] = sqrt(5.0 / 16.0) * (12.0 * y[1] * z[2] - 6.0 * x[2] * y[1] - 2.0 * y[3]);
              dYdy[2] = sqrt(5.0 / 16.0) * (12.0 * x[1] * z[2] - 2.0 * x[3] - 6.0 * x[1] * y[2]);
              dYdz[2] = sqrt(5.0 / 16.0) * (24.0 * x[1] * y[1] * z[1]);
              dYdx[3] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * y[1] * z[1]);
              dYdy[3] = sqrt(5.0 / 8.0) * (4.0 * z[3] - 3.0 * x[2] * z[1] - 9.0 * y[2] * z[1]);
              dYdz[3] = sqrt(5.0 / 8.0) * (12.0 * y[1] * z[2] - 3.0 * x[2] * y[1] - 3.0 * y[3]);
              dYdx[4] = sqrt(1.0 / 64.0) * (-48.0 * x[1] * z[2] + 12.0 * x[3] + 12.0 * x[1] * y[2]);
              dYdy[4] = sqrt(1.0 / 64.0) * (-48.0 * y[1] * z[2] + 12.0 * y[3] + 12.0 * x[2] * y[1]);
              dYdz[4] = sqrt(1.0 / 64.0) * (32.0 * z[3] - 48.0 * x[2] * z[1] - 48.0 * y[2] * z[1]);
              dYdx[5] = sqrt(5.0 / 8.0) * (4.0 * z[3] - 9.0 * x[2] * z[1] - 3.0 * y[2] * z[1]);
              dYdy[5] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * y[1] * z[1]);
              dYdz[5] = sqrt(5.0 / 8.0) * (12.0 * x[1] * z[2] - 3.0 * x[3] - 3.0 * x[1] * y[2]);
              dYdx[6] = sqrt(5.0 / 16.0) * (12.0 * x[1] * z[2] - 4.0 * x[3]);
              dYdy[6] = sqrt(5.0 / 16.0) * (-12.0 * y[1] * z[2] + 4.0 * y[3]);
              dYdz[6] = sqrt(5.0 / 16.0) * (12.0 * x[2] * z[1] - 12.0 * y[2] * z[1]);
              dYdx[7] = sqrt(35.0 / 8.0) * (3.0 * x[2] * z[1] - 3.0 * y[2] * z[1]);
              dYdy[7] = sqrt(35.0 / 8.0) * (-6.0 * x[1] * y[1] * z[1]);
              dYdz[7] = sqrt(35.0 / 8.0) * (1.0 * x[3] - 3.0 * x[1] * y[2]);
              dYdx[8] = sqrt(35.0 / 64.0) * (4.0 * x[3] - 12.0 * x[1] * y[2]);
              dYdy[8] = sqrt(35.0 / 64.0) * (-12.0 * x[2] * y[1] + 4.0 * y[3]);
              dYdz[8] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = sqrt(35.0 / 64.0) * (24.0 * x[1] * y[1]);
                d2Ydxdy[0] = sqrt(35.0 / 64.0) * (12.0 * x[2] - 12.0 * y[2]);
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = sqrt(35.0 / 64.0) * (-24.0 * x[1] * y[1]);
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = sqrt(35.0 / 8.0) * (6.0 * y[1] * z[1]);
                d2Ydxdy[1] = sqrt(35.0 / 8.0) * (6.0 * x[1] * z[1]);
                d2Ydxdz[1] = sqrt(35.0 / 8.0) * (6.0 * x[1] * y[1]);
                d2Ydydy[1] = sqrt(35.0 / 8.0) * (-6.0 * y[1] * z[1]);
                d2Ydydz[1] = sqrt(35.0 / 8.0) * (3.0 * x[2] - 3.0 * y[2]);
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = sqrt(5.0 / 16.0) * (-12.0 * x[1] * y[1]);
                d2Ydxdy[2] = sqrt(5.0 / 16.0) * (12.0 * z[2] - 6.0 * x[2] - 6.0 * y[2]);
                d2Ydxdz[2] = sqrt(5.0 / 16.0) * (24.0 * y[1] * z[1]);
                d2Ydydy[2] = sqrt(5.0 / 16.0) * (-12.0 * x[1] * y[1]);
                d2Ydydz[2] = sqrt(5.0 / 16.0) * (24.0 * x[1] * z[1]);
                d2Ydzdz[2] = sqrt(5.0 / 16.0) * (24.0 * x[1] * y[1]);
                d2Ydxdx[3] = sqrt(5.0 / 8.0) * (-6.0 * y[1] * z[1]);
                d2Ydxdy[3] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * z[1]);
                d2Ydxdz[3] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * y[1]);
                d2Ydydy[3] = sqrt(5.0 / 8.0) * (-18.0 * y[1] * z[1]);
                d2Ydydz[3] = sqrt(5.0 / 8.0) * (12.0 * z[2] - 3.0 * x[2] - 9.0 * y[2]);
                d2Ydzdz[3] = sqrt(5.0 / 8.0) * (24.0 * y[1] * z[1]);
                d2Ydxdx[4] = sqrt(1.0 / 64.0) * (-48.0 * z[2] + 36.0 * x[2] + 12.0 * y[2]);
                d2Ydxdy[4] = sqrt(1.0 / 64.0) * (24.0 * x[1] * y[1]);
                d2Ydxdz[4] = sqrt(1.0 / 64.0) * (-96.0 * x[1] * z[1]);
                d2Ydydy[4] = sqrt(1.0 / 64.0) * (-48.0 * z[2] + 36.0 * y[2] + 12.0 * x[2]);
                d2Ydydz[4] = sqrt(1.0 / 64.0) * (-96.0 * y[1] * z[1]);
                d2Ydzdz[4] = sqrt(1.0 / 64.0) * (96.0 * z[2] - 48.0 * x[2] - 48.0 * y[2]);
                d2Ydxdx[5] = sqrt(5.0 / 8.0) * (-18.0 * x[1] * z[1]);
                d2Ydxdy[5] = sqrt(5.0 / 8.0) * (-6.0 * y[1] * z[1]);
                d2Ydxdz[5] = sqrt(5.0 / 8.0) * (12.0 * z[2] - 9.0 * x[2] - 3.0 * y[2]);
                d2Ydydy[5] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * z[1]);
                d2Ydydz[5] = sqrt(5.0 / 8.0) * (-6.0 * x[1] * y[1]);
                d2Ydzdz[5] = sqrt(5.0 / 8.0) * (24.0 * x[1] * z[1]);
                d2Ydxdx[6] = sqrt(5.0 / 16.0) * (12.0 * z[2] - 12.0 * x[2]);
                d2Ydxdy[6] = 0.0;
                d2Ydxdz[6] = sqrt(5.0 / 16.0) * (24.0 * x[1] * z[1]);
                d2Ydydy[6] = sqrt(5.0 / 16.0) * (-12.0 * z[2] + 12.0 * y[2]);
                d2Ydydz[6] = sqrt(5.0 / 16.0) * (-24.0 * y[1] * z[1]);
                d2Ydzdz[6] = sqrt(5.0 / 16.0) * (12.0 * x[2] - 12.0 * y[2]);
                d2Ydxdx[7] = sqrt(35.0 / 8.0) * (6.0 * x[1] * z[1]);
                d2Ydxdy[7] = sqrt(35.0 / 8.0) * (-6.0 * y[1] * z[1]);
                d2Ydxdz[7] = sqrt(35.0 / 8.0) * (3.0 * x[2] - 3.0 * y[2]);
                d2Ydydy[7] = sqrt(35.0 / 8.0) * (-6.0 * x[1] * z[1]);
                d2Ydydz[7] = sqrt(35.0 / 8.0) * (-6.0 * x[1] * y[1]);
                d2Ydzdz[7] = 0.0;
                d2Ydxdx[8] = sqrt(35.0 / 64.0) * (12.0 * x[2] - 12.0 * y[2]);
                d2Ydxdy[8] = sqrt(35.0 / 64.0) * (-24.0 * x[1] * y[1]);
                d2Ydxdz[8] = 0.0;
                d2Ydydy[8] = sqrt(35.0 / 64.0) * (-12.0 * x[2] + 12.0 * y[2]);
                d2Ydydz[8] = 0.0;
                d2Ydzdz[8] = 0.0;
              }
            }
            break;
          case 5:
            Y[0] = sqrt(63.0 / 128.0) * (5.0 * x[4] * y[1] - 10.0 * x[2] * y[3] + 1.0 * y[5]);
            Y[1] = sqrt(315.0 / 64.0) * (4.0 * x[3] * y[1] * z[1] - 4.0 * x[1] * y[3] * z[1]);
            Y[2] = sqrt(35.0 / 128.0) *
                   (24.0 * x[2] * y[1] * z[2] - 3.0 * x[4] * y[1] - 2.0 * x[2] * y[3] - 8.0 * y[3] * z[2] + 1.0 * y[5]);
            Y[3] = sqrt(105.0 / 16.0) * (4.0 * x[1] * y[1] * z[3] - 2.0 * x[3] * y[1] * z[1] - 2.0 * x[1] * y[3] * z[1]);
            Y[4] = sqrt(15.0 / 64.0) * (8.0 * y[1] * z[4] - 12.0 * x[2] * y[1] * z[2] - 12.0 * y[3] * z[2] +
                                        1.0 * x[4] * y[1] + 1.0 * y[5] + 2.0 * x[2] * y[3]);
            Y[5] = sqrt(1.0 / 64.0) * (8.0 * z[5] - 40.0 * x[2] * z[3] - 40.0 * y[2] * z[3] + 15.0 * x[4] * z[1] +
                                       15.0 * y[4] * z[1] + 30.0 * x[2] * y[2] * z[1]);
            Y[6] = sqrt(15.0 / 64.0) * (8.0 * x[1] * z[4] - 12.0 * x[3] * z[2] - 12.0 * x[1] * y[2] * z[2] +
                                        1.0 * x[5] + 1.0 * x[1] * y[4] + 2.0 * x[3] * y[2]);
            Y[7] = sqrt(105.0 / 16.0) * (2.0 * x[2] * z[3] - 1.0 * x[4] * z[1] - 2.0 * y[2] * z[3] + 1.0 * y[4] * z[1]);
            Y[8] = sqrt(35.0 / 128.0) *
                   (8.0 * x[3] * z[2] - 1.0 * x[5] + 2.0 * x[3] * y[2] - 24.0 * x[1] * y[2] * z[2] + 3.0 * x[1] * y[4]);
            Y[9] = sqrt(315.0 / 64.0) * (1.0 * x[4] * z[1] - 6.0 * x[2] * y[2] * z[1] + 1.0 * y[4] * z[1]);
            Y[10] = sqrt(63.0 / 128.0) * (1.0 * x[5] - 10.0 * x[3] * y[2] + 5.0 * x[1] * y[4]);
            if (_highestDerivative >= 1) {
              dYdx[0] = sqrt(63.0 / 128.0) * (20.0 * x[3] * y[1] - 20.0 * x[1] * y[3]);
              dYdy[0] = sqrt(63.0 / 128.0) * (5.0 * x[4] - 30.0 * x[2] * y[2] + 5.0 * y[4]);
              dYdz[0] = 0.0;
              dYdx[1] = sqrt(315.0 / 64.0) * (12.0 * x[2] * y[1] * z[1] - 4.0 * y[3] * z[1]);
              dYdy[1] = sqrt(315.0 / 64.0) * (4.0 * x[3] * z[1] - 12.0 * x[1] * y[2] * z[1]);
              dYdz[1] = sqrt(315.0 / 64.0) * (4.0 * x[3] * y[1] - 4.0 * x[1] * y[3]);
              dYdx[2] = sqrt(35.0 / 128.0) * (48.0 * x[1] * y[1] * z[2] - 12.0 * x[3] * y[1] - 4.0 * x[1] * y[3]);
              dYdy[2] = sqrt(35.0 / 128.0) *
                        (24.0 * x[2] * z[2] - 3.0 * x[4] - 6.0 * x[2] * y[2] - 24.0 * y[2] * z[2] + 5.0 * y[4]);
              dYdz[2] = sqrt(35.0 / 128.0) * (48.0 * x[2] * y[1] * z[1] - 16.0 * y[3] * z[1]);
              dYdx[3] = sqrt(105.0 / 16.0) * (4.0 * y[1] * z[3] - 6.0 * x[2] * y[1] * z[1] - 2.0 * y[3] * z[1]);
              dYdy[3] = sqrt(105.0 / 16.0) * (4.0 * x[1] * z[3] - 2.0 * x[3] * z[1] - 6.0 * x[1] * y[2] * z[1]);
              dYdz[3] = sqrt(105.0 / 16.0) * (12.0 * x[1] * y[1] * z[2] - 2.0 * x[3] * y[1] - 2.0 * x[1] * y[3]);
              dYdx[4] = sqrt(15.0 / 64.0) * (-24.0 * x[1] * y[1] * z[2] + 4.0 * x[3] * y[1] + 4.0 * x[1] * y[3]);
              dYdy[4] = sqrt(15.0 / 64.0) * (8.0 * z[4] - 12.0 * x[2] * z[2] - 36.0 * y[2] * z[2] + 1.0 * x[4] +
                                             5.0 * y[4] + 6.0 * x[2] * y[2]);
              dYdz[4] = sqrt(15.0 / 64.0) * (32.0 * y[1] * z[3] - 24.0 * x[2] * y[1] * z[1] - 24.0 * y[3] * z[1]);
              dYdx[5] = sqrt(1.0 / 64.0) * (-80.0 * x[1] * z[3] + 60.0 * x[3] * z[1] + 60.0 * x[1] * y[2] * z[1]);
              dYdy[5] = sqrt(1.0 / 64.0) * (-80.0 * y[1] * z[3] + 60.0 * y[3] * z[1] + 60.0 * x[2] * y[1] * z[1]);
              dYdz[5] = sqrt(1.0 / 64.0) * (40.0 * z[4] - 120.0 * x[2] * z[2] - 120.0 * y[2] * z[2] + 15.0 * x[4] +
                                            15.0 * y[4] + 30.0 * x[2] * y[2]);
              dYdx[6] = sqrt(15.0 / 64.0) * (8.0 * z[4] - 36.0 * x[2] * z[2] - 12.0 * y[2] * z[2] + 5.0 * x[4] +
                                             1.0 * y[4] + 6.0 * x[2] * y[2]);
              dYdy[6] = sqrt(15.0 / 64.0) * (-24.0 * x[1] * y[1] * z[2] + 4.0 * x[1] * y[3] + 4.0 * x[3] * y[1]);
              dYdz[6] = sqrt(15.0 / 64.0) * (32.0 * x[1] * z[3] - 24.0 * x[3] * z[1] - 24.0 * x[1] * y[2] * z[1]);
              dYdx[7] = sqrt(105.0 / 16.0) * (4.0 * x[1] * z[3] - 4.0 * x[3] * z[1]);
              dYdy[7] = sqrt(105.0 / 16.0) * (-4.0 * y[1] * z[3] + 4.0 * y[3] * z[1]);
              dYdz[7] = sqrt(105.0 / 16.0) * (6.0 * x[2] * z[2] - 1.0 * x[4] - 6.0 * y[2] * z[2] + 1.0 * y[4]);
              dYdx[8] = sqrt(35.0 / 128.0) *
                        (24.0 * x[2] * z[2] - 5.0 * x[4] + 6.0 * x[2] * y[2] - 24.0 * y[2] * z[2] + 3.0 * y[4]);
              dYdy[8] = sqrt(35.0 / 128.0) * (4.0 * x[3] * y[1] - 48.0 * x[1] * y[1] * z[2] + 12.0 * x[1] * y[3]);
              dYdz[8] = sqrt(35.0 / 128.0) * (16.0 * x[3] * z[1] - 48.0 * x[1] * y[2] * z[1]);
              dYdx[9] = sqrt(315.0 / 64.0) * (4.0 * x[3] * z[1] - 12.0 * x[1] * y[2] * z[1]);
              dYdy[9] = sqrt(315.0 / 64.0) * (-12.0 * x[2] * y[1] * z[1] + 4.0 * y[3] * z[1]);
              dYdz[9] = sqrt(315.0 / 64.0) * (1.0 * x[4] - 6.0 * x[2] * y[2] + 1.0 * y[4]);
              dYdx[10] = sqrt(63.0 / 128.0) * (5.0 * x[4] - 30.0 * x[2] * y[2] + 5.0 * y[4]);
              dYdy[10] = sqrt(63.0 / 128.0) * (-20.0 * x[3] * y[1] + 20.0 * x[1] * y[3]);
              dYdz[10] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = sqrt(63.0 / 128.0) * (60.0 * x[2] * y[1] - 20.0 * y[3]);
                d2Ydxdy[0] = sqrt(63.0 / 128.0) * (20.0 * x[3] - 60.0 * x[1] * y[2]);
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = sqrt(63.0 / 128.0) * (-60.0 * x[2] * y[1] + 20.0 * y[3]);
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = sqrt(315.0 / 64.0) * (24.0 * x[1] * y[1] * z[1]);
                d2Ydxdy[1] = sqrt(315.0 / 64.0) * (12.0 * x[2] * z[1] - 12.0 * y[2] * z[1]);
                d2Ydxdz[1] = sqrt(315.0 / 64.0) * (12.0 * x[2] * y[1] - 4.0 * y[3]);
                d2Ydydy[1] = sqrt(315.0 / 64.0) * (-24.0 * x[1] * y[1] * z[1]);
                d2Ydydz[1] = sqrt(315.0 / 64.0) * (4.0 * x[3] - 12.0 * x[1] * y[2]);
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = sqrt(35.0 / 128.0) * (48.0 * y[1] * z[2] - 36.0 * x[2] * y[1] - 4.0 * y[3]);
                d2Ydxdy[2] = sqrt(35.0 / 128.0) * (48.0 * x[1] * z[2] - 12.0 * x[3] - 12.0 * x[1] * y[2]);
                d2Ydxdz[2] = sqrt(35.0 / 128.0) * (96.0 * x[1] * y[1] * z[1]);
                d2Ydydy[2] = sqrt(35.0 / 128.0) * (-12.0 * x[2] * y[1] - 48.0 * y[1] * z[2] + 20.0 * y[3]);
                d2Ydydz[2] = sqrt(35.0 / 128.0) * (48.0 * x[2] * z[1] - 48.0 * y[2] * z[1]);
                d2Ydzdz[2] = sqrt(35.0 / 128.0) * (48.0 * x[2] * y[1] - 16.0 * y[3]);
                d2Ydxdx[3] = sqrt(105.0 / 16.0) * (-12.0 * x[1] * y[1] * z[1]);
                d2Ydxdy[3] = sqrt(105.0 / 16.0) * (4.0 * z[3] - 6.0 * x[2] * z[1] - 6.0 * y[2] * z[1]);
                d2Ydxdz[3] = sqrt(105.0 / 16.0) * (12.0 * y[1] * z[2] - 6.0 * x[2] * y[1] - 2.0 * y[3]);
                d2Ydydy[3] = sqrt(105.0 / 16.0) * (-12.0 * x[1] * y[1] * z[1]);
                d2Ydydz[3] = sqrt(105.0 / 16.0) * (12.0 * x[1] * z[2] - 2.0 * x[3] - 6.0 * x[1] * y[2]);
                d2Ydzdz[3] = sqrt(105.0 / 16.0) * (24.0 * x[1] * y[1] * z[1]);
                d2Ydxdx[4] = sqrt(15.0 / 64.0) * (-24.0 * y[1] * z[2] + 12.0 * x[2] * y[1] + 4.0 * y[3]);
                d2Ydxdy[4] = sqrt(15.0 / 64.0) * (-24.0 * x[1] * z[2] + 4.0 * x[3] + 12.0 * x[1] * y[2]);
                d2Ydxdz[4] = sqrt(15.0 / 64.0) * (-48.0 * x[1] * y[1] * z[1]);
                d2Ydydy[4] = sqrt(15.0 / 64.0) * (-72.0 * y[1] * z[2] + 20.0 * y[3] + 12.0 * x[2] * y[1]);
                d2Ydydz[4] = sqrt(15.0 / 64.0) * (32.0 * z[3] - 24.0 * x[2] * z[1] - 72.0 * y[2] * z[1]);
                d2Ydzdz[4] = sqrt(15.0 / 64.0) * (96.0 * y[1] * z[2] - 24.0 * x[2] * y[1] - 24.0 * y[3]);
                d2Ydxdx[5] = sqrt(1.0 / 64.0) * (-80.0 * z[3] + 180.0 * x[2] * z[1] + 60.0 * y[2] * z[1]);
                d2Ydxdy[5] = sqrt(1.0 / 64.0) * (120.0 * x[1] * y[1] * z[1]);
                d2Ydxdz[5] = sqrt(1.0 / 64.0) * (-240.0 * x[1] * z[2] + 60.0 * x[3] + 60.0 * x[1] * y[2]);
                d2Ydydy[5] = sqrt(1.0 / 64.0) * (-80.0 * z[3] + 180.0 * y[2] * z[1] + 60.0 * x[2] * z[1]);
                d2Ydydz[5] = sqrt(1.0 / 64.0) * (-240.0 * y[1] * z[2] + 60.0 * y[3] + 60.0 * x[2] * y[1]);
                d2Ydzdz[5] = sqrt(1.0 / 64.0) * (160.0 * z[3] - 240.0 * x[2] * z[1] - 240.0 * y[2] * z[1]);
                d2Ydxdx[6] = sqrt(15.0 / 64.0) * (-72.0 * x[1] * z[2] + 20.0 * x[3] + 12.0 * x[1] * y[2]);
                d2Ydxdy[6] = sqrt(15.0 / 64.0) * (-24.0 * y[1] * z[2] + 4.0 * y[3] + 12.0 * x[2] * y[1]);
                d2Ydxdz[6] = sqrt(15.0 / 64.0) * (32.0 * z[3] - 72.0 * x[2] * z[1] - 24.0 * y[2] * z[1]);
                d2Ydydy[6] = sqrt(15.0 / 64.0) * (-24.0 * x[1] * z[2] + 12.0 * x[1] * y[2] + 4.0 * x[3]);
                d2Ydydz[6] = sqrt(15.0 / 64.0) * (-48.0 * x[1] * y[1] * z[1]);
                d2Ydzdz[6] = sqrt(15.0 / 64.0) * (96.0 * x[1] * z[2] - 24.0 * x[3] - 24.0 * x[1] * y[2]);
                d2Ydxdx[7] = sqrt(105.0 / 16.0) * (4.0 * z[3] - 12.0 * x[2] * z[1]);
                d2Ydxdy[7] = 0.0;
                d2Ydxdz[7] = sqrt(105.0 / 16.0) * (12.0 * x[1] * z[2] - 4.0 * x[3]);
                d2Ydydy[7] = sqrt(105.0 / 16.0) * (-4.0 * z[3] + 12.0 * y[2] * z[1]);
                d2Ydydz[7] = sqrt(105.0 / 16.0) * (-12.0 * y[1] * z[2] + 4.0 * y[3]);
                d2Ydzdz[7] = sqrt(105.0 / 16.0) * (12.0 * x[2] * z[1] - 12.0 * y[2] * z[1]);
                d2Ydxdx[8] = sqrt(35.0 / 128.0) * (48.0 * x[1] * z[2] - 20.0 * x[3] + 12.0 * x[1] * y[2]);
                d2Ydxdy[8] = sqrt(35.0 / 128.0) * (12.0 * x[2] * y[1] - 48.0 * y[1] * z[2] + 12.0 * y[3]);
                d2Ydxdz[8] = sqrt(35.0 / 128.0) * (48.0 * x[2] * z[1] - 48.0 * y[2] * z[1]);
                d2Ydydy[8] = sqrt(35.0 / 128.0) * (4.0 * x[3] - 48.0 * x[1] * z[2] + 36.0 * x[1] * y[2]);
                d2Ydydz[8] = sqrt(35.0 / 128.0) * (-96.0 * x[1] * y[1] * z[1]);
                d2Ydzdz[8] = sqrt(35.0 / 128.0) * (16.0 * x[3] - 48.0 * x[1] * y[2]);
                d2Ydxdx[9] = sqrt(315.0 / 64.0) * (12.0 * x[2] * z[1] - 12.0 * y[2] * z[1]);
                d2Ydxdy[9] = sqrt(315.0 / 64.0) * (-24.0 * x[1] * y[1] * z[1]);
                d2Ydxdz[9] = sqrt(315.0 / 64.0) * (4.0 * x[3] - 12.0 * x[1] * y[2]);
                d2Ydydy[9] = sqrt(315.0 / 64.0) * (-12.0 * x[2] * z[1] + 12.0 * y[2] * z[1]);
                d2Ydydz[9] = sqrt(315.0 / 64.0) * (-12.0 * x[2] * y[1] + 4.0 * y[3]);
                d2Ydzdz[9] = 0.0;
                d2Ydxdx[10] = sqrt(63.0 / 128.0) * (20.0 * x[3] - 60.0 * x[1] * y[2]);
                d2Ydxdy[10] = sqrt(63.0 / 128.0) * (-60.0 * x[2] * y[1] + 20.0 * y[3]);
                d2Ydxdz[10] = 0.0;
                d2Ydydy[10] = sqrt(63.0 / 128.0) * (-20.0 * x[3] + 60.0 * x[1] * y[2]);
                d2Ydydz[10] = 0.0;
                d2Ydzdz[10] = 0.0;
              }
            }
            break;
          case 6:
            Y[0] = sqrt(231.0 / 512.0) * (6.0 * x[5] * y[1] - 20.0 * x[3] * y[3] + 6.0 * x[1] * y[5]);
            Y[1] = sqrt(693.0 / 128.0) * (5.0 * x[4] * y[1] * z[1] - 10.0 * x[2] * y[3] * z[1] + 1.0 * y[5] * z[1]);
            Y[2] = sqrt(63.0 / 256.0) *
                   (40.0 * x[3] * y[1] * z[2] - 4.0 * x[5] * y[1] - 40.0 * x[1] * y[3] * z[2] + 4.0 * x[1] * y[5]);
            Y[3] = sqrt(105.0 / 128.0) * (24.0 * x[2] * y[1] * z[3] - 9.0 * x[4] * y[1] * z[1] -
                                          6.0 * x[2] * y[3] * z[1] - 8.0 * y[3] * z[3] + 3.0 * y[5] * z[1]);
            Y[4] = sqrt(105.0 / 512.0) * (32.0 * x[1] * y[1] * z[4] - 32.0 * x[3] * y[1] * z[2] - 32.0 * x[1] * y[3] * z[2] +
                                          2.0 * x[5] * y[1] + 2.0 * x[1] * y[5] + 4.0 * x[3] * y[3]);
            Y[5] = sqrt(21.0 / 64.0) * (8.0 * y[1] * z[5] - 20.0 * x[2] * y[1] * z[3] - 20.0 * y[3] * z[3] +
                                        5.0 * x[4] * y[1] * z[1] + 5.0 * y[5] * z[1] + 10.0 * x[2] * y[3] * z[1]);
            Y[6] = sqrt(1.0 / 256.0) * (331.0 * z[6] - 120.0 * x[2] * z[4] - 120.0 * y[2] * z[4] - 315.0 * z[6] +
                                        90.0 * x[4] * z[2] + 90.0 * y[4] * z[2] + 180.0 * x[2] * y[2] * z[2] -
                                        5.0 * x[6] - 5.0 * y[6] - 15.0 * x[4] * y[2] - 15.0 * x[2] * y[4]);
            Y[7] = sqrt(21.0 / 64.0) * (8.0 * x[1] * z[5] - 20.0 * x[3] * z[3] - 20.0 * x[1] * y[2] * z[3] +
                                        5.0 * x[5] * z[1] + 5.0 * x[1] * y[4] * z[1] + 10.0 * x[3] * y[2] * z[1]);
            Y[8] = sqrt(105.0 / 512.0) *
                   (16.0 * x[2] * z[4] - 16.0 * x[4] * z[2] + 1.0 * x[6] - 1.0 * x[2] * y[4] + 1.0 * x[4] * y[2] -
                    34.0 * y[2] * z[4] + 16.0 * y[4] * z[2] + 18.0 * y[2] * z[4] - 1.0 * y[6]);
            Y[9] = sqrt(105.0 / 128.0) * (8.0 * x[3] * z[3] - 3.0 * x[5] * z[1] + 6.0 * x[3] * y[2] * z[1] -
                                          24.0 * x[1] * y[2] * z[3] + 9.0 * x[1] * y[4] * z[1]);
            Y[10] = sqrt(63.0 / 256.0) * (10.0 * x[4] * z[2] - 1.0 * x[6] + 5.0 * x[4] * y[2] - 60.0 * x[2] * y[2] * z[2] +
                                          5.0 * x[2] * y[4] + 10.0 * y[4] * z[2] - 1.0 * y[6]);
            Y[11] = sqrt(693.0 / 128.0) * (1.0 * x[5] * z[1] - 10.0 * x[3] * y[2] * z[1] + 5.0 * x[1] * y[4] * z[1]);
            Y[12] = sqrt(231.0 / 512.0) * (1.0 * x[6] - 15.0 * x[4] * y[2] + 15.0 * x[2] * y[4] - 1.0 * y[6]);
            if (_highestDerivative >= 1) {
              dYdx[0] = sqrt(231.0 / 512.0) * (30.0 * x[4] * y[1] - 60.0 * x[2] * y[3] + 6.0 * y[5]);
              dYdy[0] = sqrt(231.0 / 512.0) * (6.0 * x[5] - 60.0 * x[3] * y[2] + 30.0 * x[1] * y[4]);
              dYdz[0] = 0.0;
              dYdx[1] = sqrt(693.0 / 128.0) * (20.0 * x[3] * y[1] * z[1] - 20.0 * x[1] * y[3] * z[1]);
              dYdy[1] = sqrt(693.0 / 128.0) * (5.0 * x[4] * z[1] - 30.0 * x[2] * y[2] * z[1] + 5.0 * y[4] * z[1]);
              dYdz[1] = sqrt(693.0 / 128.0) * (5.0 * x[4] * y[1] - 10.0 * x[2] * y[3] + 1.0 * y[5]);
              dYdx[2] =
                  sqrt(63.0 / 256.0) * (120.0 * x[2] * y[1] * z[2] - 20.0 * x[4] * y[1] - 40.0 * y[3] * z[2] + 4.0 * y[5]);
              dYdy[2] =
                  sqrt(63.0 / 256.0) * (40.0 * x[3] * z[2] - 4.0 * x[5] - 120.0 * x[1] * y[2] * z[2] + 20.0 * x[1] * y[4]);
              dYdz[2] = sqrt(63.0 / 256.0) * (80.0 * x[3] * y[1] * z[1] - 80.0 * x[1] * y[3] * z[1]);
              dYdx[3] = sqrt(105.0 / 128.0) *
                        (48.0 * x[1] * y[1] * z[3] - 36.0 * x[3] * y[1] * z[1] - 12.0 * x[1] * y[3] * z[1]);
              dYdy[3] = sqrt(105.0 / 128.0) * (24.0 * x[2] * z[3] - 9.0 * x[4] * z[1] - 18.0 * x[2] * y[2] * z[1] -
                                               24.0 * y[2] * z[3] + 15.0 * y[4] * z[1]);
              dYdz[3] = sqrt(105.0 / 128.0) * (72.0 * x[2] * y[1] * z[2] - 9.0 * x[4] * y[1] - 6.0 * x[2] * y[3] -
                                               24.0 * y[3] * z[2] + 3.0 * y[5]);
              dYdx[4] = sqrt(105.0 / 512.0) * (32.0 * y[1] * z[4] - 96.0 * x[2] * y[1] * z[2] - 32.0 * y[3] * z[2] +
                                               10.0 * x[4] * y[1] + 2.0 * y[5] + 12.0 * x[2] * y[3]);
              dYdy[4] = sqrt(105.0 / 512.0) * (32.0 * x[1] * z[4] - 32.0 * x[3] * z[2] - 96.0 * x[1] * y[2] * z[2] +
                                               2.0 * x[5] + 10.0 * x[1] * y[4] + 12.0 * x[3] * y[2]);
              dYdz[4] = sqrt(105.0 / 512.0) *
                        (128.0 * x[1] * y[1] * z[3] - 64.0 * x[3] * y[1] * z[1] - 64.0 * x[1] * y[3] * z[1]);
              dYdx[5] =
                  sqrt(21.0 / 64.0) * (-40.0 * x[1] * y[1] * z[3] + 20.0 * x[3] * y[1] * z[1] + 20.0 * x[1] * y[3] * z[1]);
              dYdy[5] = sqrt(21.0 / 64.0) * (8.0 * z[5] - 20.0 * x[2] * z[3] - 60.0 * y[2] * z[3] + 5.0 * x[4] * z[1] +
                                             25.0 * y[4] * z[1] + 30.0 * x[2] * y[2] * z[1]);
              dYdz[5] = sqrt(21.0 / 64.0) * (40.0 * y[1] * z[4] - 60.0 * x[2] * y[1] * z[2] - 60.0 * y[3] * z[2] +
                                             5.0 * x[4] * y[1] + 5.0 * y[5] + 10.0 * x[2] * y[3]);
              dYdx[6] = sqrt(1.0 / 256.0) * (-240.0 * x[1] * z[4] + 360.0 * x[3] * z[2] + 360.0 * x[1] * y[2] * z[2] -
                                             30.0 * x[5] - 60.0 * x[3] * y[2] - 30.0 * x[1] * y[4]);
              dYdy[6] = sqrt(1.0 / 256.0) * (-240.0 * y[1] * z[4] + 360.0 * y[3] * z[2] + 360.0 * x[2] * y[1] * z[2] -
                                             30.0 * y[5] - 30.0 * x[4] * y[1] - 60.0 * x[2] * y[3]);
              dYdz[6] = sqrt(1.0 / 256.0) * (1986.0 * z[5] - 480.0 * x[2] * z[3] - 480.0 * y[2] * z[3] - 1890.0 * z[5] +
                                             180.0 * x[4] * z[1] + 180.0 * y[4] * z[1] + 360.0 * x[2] * y[2] * z[1]);
              dYdx[7] = sqrt(21.0 / 64.0) * (8.0 * z[5] - 60.0 * x[2] * z[3] - 20.0 * y[2] * z[3] + 25.0 * x[4] * z[1] +
                                             5.0 * y[4] * z[1] + 30.0 * x[2] * y[2] * z[1]);
              dYdy[7] =
                  sqrt(21.0 / 64.0) * (-40.0 * x[1] * y[1] * z[3] + 20.0 * x[1] * y[3] * z[1] + 20.0 * x[3] * y[1] * z[1]);
              dYdz[7] = sqrt(21.0 / 64.0) * (40.0 * x[1] * z[4] - 60.0 * x[3] * z[2] - 60.0 * x[1] * y[2] * z[2] +
                                             5.0 * x[5] + 5.0 * x[1] * y[4] + 10.0 * x[3] * y[2]);
              dYdx[8] = sqrt(105.0 / 512.0) *
                        (32.0 * x[1] * z[4] - 64.0 * x[3] * z[2] + 6.0 * x[5] - 2.0 * x[1] * y[4] + 4.0 * x[3] * y[2]);
              dYdy[8] = sqrt(105.0 / 512.0) * (-4.0 * x[2] * y[3] + 2.0 * x[4] * y[1] - 68.0 * y[1] * z[4] +
                                               64.0 * y[3] * z[2] + 36.0 * y[1] * z[4] - 6.0 * y[5]);
              dYdz[8] = sqrt(105.0 / 512.0) * (64.0 * x[2] * z[3] - 32.0 * x[4] * z[1] - 136.0 * y[2] * z[3] +
                                               32.0 * y[4] * z[1] + 72.0 * y[2] * z[3]);
              dYdx[9] = sqrt(105.0 / 128.0) * (24.0 * x[2] * z[3] - 15.0 * x[4] * z[1] + 18.0 * x[2] * y[2] * z[1] -
                                               24.0 * y[2] * z[3] + 9.0 * y[4] * z[1]);
              dYdy[9] = sqrt(105.0 / 128.0) *
                        (12.0 * x[3] * y[1] * z[1] - 48.0 * x[1] * y[1] * z[3] + 36.0 * x[1] * y[3] * z[1]);
              dYdz[9] = sqrt(105.0 / 128.0) * (24.0 * x[3] * z[2] - 3.0 * x[5] + 6.0 * x[3] * y[2] -
                                               72.0 * x[1] * y[2] * z[2] + 9.0 * x[1] * y[4]);
              dYdx[10] = sqrt(63.0 / 256.0) * (40.0 * x[3] * z[2] - 6.0 * x[5] + 20.0 * x[3] * y[2] -
                                               120.0 * x[1] * y[2] * z[2] + 10.0 * x[1] * y[4]);
              dYdy[10] = sqrt(63.0 / 256.0) * (10.0 * x[4] * y[1] - 120.0 * x[2] * y[1] * z[2] + 20.0 * x[2] * y[3] +
                                               40.0 * y[3] * z[2] - 6.0 * y[5]);
              dYdz[10] = sqrt(63.0 / 256.0) * (20.0 * x[4] * z[1] - 120.0 * x[2] * y[2] * z[1] + 20.0 * y[4] * z[1]);
              dYdx[11] = sqrt(693.0 / 128.0) * (5.0 * x[4] * z[1] - 30.0 * x[2] * y[2] * z[1] + 5.0 * y[4] * z[1]);
              dYdy[11] = sqrt(693.0 / 128.0) * (-20.0 * x[3] * y[1] * z[1] + 20.0 * x[1] * y[3] * z[1]);
              dYdz[11] = sqrt(693.0 / 128.0) * (1.0 * x[5] - 10.0 * x[3] * y[2] + 5.0 * x[1] * y[4]);
              dYdx[12] = sqrt(231.0 / 512.0) * (6.0 * x[5] - 60.0 * x[3] * y[2] + 30.0 * x[1] * y[4]);
              dYdy[12] = sqrt(231.0 / 512.0) * (-30.0 * x[4] * y[1] + 60.0 * x[2] * y[3] - 6.0 * y[5]);
              dYdz[12] = 0.0;
              if (_highestDerivative >= 2) {
                d2Ydxdx[0] = sqrt(231.0 / 512.0) * (120.0 * x[3] * y[1] - 120.0 * x[1] * y[3]);
                d2Ydxdy[0] = sqrt(231.0 / 512.0) * (30.0 * x[4] - 180.0 * x[2] * y[2] + 30.0 * y[4]);
                d2Ydxdz[0] = 0.0;
                d2Ydydy[0] = sqrt(231.0 / 512.0) * (-120.0 * x[3] * y[1] + 120.0 * x[1] * y[3]);
                d2Ydydz[0] = 0.0;
                d2Ydzdz[0] = 0.0;
                d2Ydxdx[1] = sqrt(693.0 / 128.0) * (60.0 * x[2] * y[1] * z[1] - 20.0 * y[3] * z[1]);
                d2Ydxdy[1] = sqrt(693.0 / 128.0) * (20.0 * x[3] * z[1] - 60.0 * x[1] * y[2] * z[1]);
                d2Ydxdz[1] = sqrt(693.0 / 128.0) * (20.0 * x[3] * y[1] - 20.0 * x[1] * y[3]);
                d2Ydydy[1] = sqrt(693.0 / 128.0) * (-60.0 * x[2] * y[1] * z[1] + 20.0 * y[3] * z[1]);
                d2Ydydz[1] = sqrt(693.0 / 128.0) * (5.0 * x[4] - 30.0 * x[2] * y[2] + 5.0 * y[4]);
                d2Ydzdz[1] = 0.0;
                d2Ydxdx[2] = sqrt(63.0 / 256.0) * (240.0 * x[1] * y[1] * z[2] - 80.0 * x[3] * y[1]);
                d2Ydxdy[2] = sqrt(63.0 / 256.0) * (120.0 * x[2] * z[2] - 20.0 * x[4] - 120.0 * y[2] * z[2] + 20.0 * y[4]);
                d2Ydxdz[2] = sqrt(63.0 / 256.0) * (240.0 * x[2] * y[1] * z[1] - 80.0 * y[3] * z[1]);
                d2Ydydy[2] = sqrt(63.0 / 256.0) * (-240.0 * x[1] * y[1] * z[2] + 80.0 * x[1] * y[3]);
                d2Ydydz[2] = sqrt(63.0 / 256.0) * (80.0 * x[3] * z[1] - 240.0 * x[1] * y[2] * z[1]);
                d2Ydzdz[2] = sqrt(63.0 / 256.0) * (80.0 * x[3] * y[1] - 80.0 * x[1] * y[3]);
                d2Ydxdx[3] = sqrt(105.0 / 128.0) * (48.0 * y[1] * z[3] - 108.0 * x[2] * y[1] * z[1] - 12.0 * y[3] * z[1]);
                d2Ydxdy[3] = sqrt(105.0 / 128.0) * (48.0 * x[1] * z[3] - 36.0 * x[3] * z[1] - 36.0 * x[1] * y[2] * z[1]);
                d2Ydxdz[3] = sqrt(105.0 / 128.0) * (144.0 * x[1] * y[1] * z[2] - 36.0 * x[3] * y[1] - 12.0 * x[1] * y[3]);
                d2Ydydy[3] = sqrt(105.0 / 128.0) * (-36.0 * x[2] * y[1] * z[1] - 48.0 * y[1] * z[3] + 60.0 * y[3] * z[1]);
                d2Ydydz[3] = sqrt(105.0 / 128.0) *
                             (72.0 * x[2] * z[2] - 9.0 * x[4] - 18.0 * x[2] * y[2] - 72.0 * y[2] * z[2] + 15.0 * y[4]);
                d2Ydzdz[3] = sqrt(105.0 / 128.0) * (144.0 * x[2] * y[1] * z[1] - 48.0 * y[3] * z[1]);
                d2Ydxdx[4] = sqrt(105.0 / 512.0) * (-192.0 * x[1] * y[1] * z[2] + 40.0 * x[3] * y[1] + 24.0 * x[1] * y[3]);
                d2Ydxdy[4] = sqrt(105.0 / 512.0) * (32.0 * z[4] - 96.0 * x[2] * z[2] - 96.0 * y[2] * z[2] +
                                                    10.0 * x[4] + 10.0 * y[4] + 36.0 * x[2] * y[2]);
                d2Ydxdz[4] = sqrt(105.0 / 512.0) * (128.0 * y[1] * z[3] - 192.0 * x[2] * y[1] * z[1] - 64.0 * y[3] * z[1]);
                d2Ydydy[4] = sqrt(105.0 / 512.0) * (-192.0 * x[1] * y[1] * z[2] + 40.0 * x[1] * y[3] + 24.0 * x[3] * y[1]);
                d2Ydydz[4] = sqrt(105.0 / 512.0) * (128.0 * x[1] * z[3] - 64.0 * x[3] * z[1] - 192.0 * x[1] * y[2] * z[1]);
                d2Ydzdz[4] = sqrt(105.0 / 512.0) * (384.0 * x[1] * y[1] * z[2] - 64.0 * x[3] * y[1] - 64.0 * x[1] * y[3]);
                d2Ydxdx[5] = sqrt(21.0 / 64.0) * (-40.0 * y[1] * z[3] + 60.0 * x[2] * y[1] * z[1] + 20.0 * y[3] * z[1]);
                d2Ydxdy[5] = sqrt(21.0 / 64.0) * (-40.0 * x[1] * z[3] + 20.0 * x[3] * z[1] + 60.0 * x[1] * y[2] * z[1]);
                d2Ydxdz[5] = sqrt(21.0 / 64.0) * (-120.0 * x[1] * y[1] * z[2] + 20.0 * x[3] * y[1] + 20.0 * x[1] * y[3]);
                d2Ydydy[5] = sqrt(21.0 / 64.0) * (-120.0 * y[1] * z[3] + 100.0 * y[3] * z[1] + 60.0 * x[2] * y[1] * z[1]);
                d2Ydydz[5] = sqrt(21.0 / 64.0) * (40.0 * z[4] - 60.0 * x[2] * z[2] - 180.0 * y[2] * z[2] + 5.0 * x[4] +
                                                  25.0 * y[4] + 30.0 * x[2] * y[2]);
                d2Ydzdz[5] = sqrt(21.0 / 64.0) * (160.0 * y[1] * z[3] - 120.0 * x[2] * y[1] * z[1] - 120.0 * y[3] * z[1]);
                d2Ydxdx[6] = sqrt(1.0 / 256.0) * (-240.0 * z[4] + 1080.0 * x[2] * z[2] + 360.0 * y[2] * z[2] -
                                                  150.0 * x[4] - 180.0 * x[2] * y[2] - 30.0 * y[4]);
                d2Ydxdy[6] = sqrt(1.0 / 256.0) * (720.0 * x[1] * y[1] * z[2] - 120.0 * x[3] * y[1] - 120.0 * x[1] * y[3]);
                d2Ydxdz[6] = sqrt(1.0 / 256.0) * (-960.0 * x[1] * z[3] + 720.0 * x[3] * z[1] + 720.0 * x[1] * y[2] * z[1]);
                d2Ydydy[6] = sqrt(1.0 / 256.0) * (-240.0 * z[4] + 1080.0 * y[2] * z[2] + 360.0 * x[2] * z[2] -
                                                  150.0 * y[4] - 30.0 * x[4] - 180.0 * x[2] * y[2]);
                d2Ydydz[6] = sqrt(1.0 / 256.0) * (-960.0 * y[1] * z[3] + 720.0 * y[3] * z[1] + 720.0 * x[2] * y[1] * z[1]);
                d2Ydzdz[6] = sqrt(1.0 / 256.0) * (9930.0 * z[4] - 1440.0 * x[2] * z[2] - 1440.0 * y[2] * z[2] -
                                                  9450.0 * z[4] + 180.0 * x[4] + 180.0 * y[4] + 360.0 * x[2] * y[2]);
                d2Ydxdx[7] = sqrt(21.0 / 64.0) * (-120.0 * x[1] * z[3] + 100.0 * x[3] * z[1] + 60.0 * x[1] * y[2] * z[1]);
                d2Ydxdy[7] = sqrt(21.0 / 64.0) * (-40.0 * y[1] * z[3] + 20.0 * y[3] * z[1] + 60.0 * x[2] * y[1] * z[1]);
                d2Ydxdz[7] = sqrt(21.0 / 64.0) * (40.0 * z[4] - 180.0 * x[2] * z[2] - 60.0 * y[2] * z[2] + 25.0 * x[4] +
                                                  5.0 * y[4] + 30.0 * x[2] * y[2]);
                d2Ydydy[7] = sqrt(21.0 / 64.0) * (-40.0 * x[1] * z[3] + 60.0 * x[1] * y[2] * z[1] + 20.0 * x[3] * z[1]);
                d2Ydydz[7] = sqrt(21.0 / 64.0) * (-120.0 * x[1] * y[1] * z[2] + 20.0 * x[1] * y[3] + 20.0 * x[3] * y[1]);
                d2Ydzdz[7] = sqrt(21.0 / 64.0) * (160.0 * x[1] * z[3] - 120.0 * x[3] * z[1] - 120.0 * x[1] * y[2] * z[1]);
                d2Ydxdx[8] = sqrt(105.0 / 512.0) *
                             (32.0 * z[4] - 192.0 * x[2] * z[2] + 30.0 * x[4] - 2.0 * y[4] + 12.0 * x[2] * y[2]);
                d2Ydxdy[8] = sqrt(105.0 / 512.0) * (-8.0 * x[1] * y[3] + 8.0 * x[3] * y[1]);
                d2Ydxdz[8] = sqrt(105.0 / 512.0) * (128.0 * x[1] * z[3] - 128.0 * x[3] * z[1]);
                d2Ydydy[8] = sqrt(105.0 / 512.0) * (-12.0 * x[2] * y[2] + 2.0 * x[4] - 68.0 * z[4] +
                                                    192.0 * y[2] * z[2] + 36.0 * z[4] - 30.0 * y[4]);
                d2Ydydz[8] = sqrt(105.0 / 512.0) * (-272.0 * y[1] * z[3] + 128.0 * y[3] * z[1] + 144.0 * y[1] * z[3]);
                d2Ydzdz[8] = sqrt(105.0 / 512.0) *
                             (192.0 * x[2] * z[2] - 32.0 * x[4] - 408.0 * y[2] * z[2] + 32.0 * y[4] + 216.0 * y[2] * z[2]);
                d2Ydxdx[9] = sqrt(105.0 / 128.0) * (48.0 * x[1] * z[3] - 60.0 * x[3] * z[1] + 36.0 * x[1] * y[2] * z[1]);
                d2Ydxdy[9] = sqrt(105.0 / 128.0) * (36.0 * x[2] * y[1] * z[1] - 48.0 * y[1] * z[3] + 36.0 * y[3] * z[1]);
                d2Ydxdz[9] = sqrt(105.0 / 128.0) *
                             (72.0 * x[2] * z[2] - 15.0 * x[4] + 18.0 * x[2] * y[2] - 72.0 * y[2] * z[2] + 9.0 * y[4]);
                d2Ydydy[9] = sqrt(105.0 / 128.0) * (12.0 * x[3] * z[1] - 48.0 * x[1] * z[3] + 108.0 * x[1] * y[2] * z[1]);
                d2Ydydz[9] = sqrt(105.0 / 128.0) * (12.0 * x[3] * y[1] - 144.0 * x[1] * y[1] * z[2] + 36.0 * x[1] * y[3]);
                d2Ydzdz[9] = sqrt(105.0 / 128.0) * (48.0 * x[3] * z[1] - 144.0 * x[1] * y[2] * z[1]);
                d2Ydxdx[10] = sqrt(63.0 / 256.0) *
                              (120.0 * x[2] * z[2] - 30.0 * x[4] + 60.0 * x[2] * y[2] - 120.0 * y[2] * z[2] + 10.0 * y[4]);
                d2Ydxdy[10] = sqrt(63.0 / 256.0) * (40.0 * x[3] * y[1] - 240.0 * x[1] * y[1] * z[2] + 40.0 * x[1] * y[3]);
                d2Ydxdz[10] = sqrt(63.0 / 256.0) * (80.0 * x[3] * z[1] - 240.0 * x[1] * y[2] * z[1]);
                d2Ydydy[10] = sqrt(63.0 / 256.0) *
                              (10.0 * x[4] - 120.0 * x[2] * z[2] + 60.0 * x[2] * y[2] + 120.0 * y[2] * z[2] - 30.0 * y[4]);
                d2Ydydz[10] = sqrt(63.0 / 256.0) * (-240.0 * x[2] * y[1] * z[1] + 80.0 * y[3] * z[1]);
                d2Ydzdz[10] = sqrt(63.0 / 256.0) * (20.0 * x[4] - 120.0 * x[2] * y[2] + 20.0 * y[4]);
                d2Ydxdx[11] = sqrt(693.0 / 128.0) * (20.0 * x[3] * z[1] - 60.0 * x[1] * y[2] * z[1]);
                d2Ydxdy[11] = sqrt(693.0 / 128.0) * (-60.0 * x[2] * y[1] * z[1] + 20.0 * y[3] * z[1]);
                d2Ydxdz[11] = sqrt(693.0 / 128.0) * (5.0 * x[4] - 30.0 * x[2] * y[2] + 5.0 * y[4]);
                d2Ydydy[11] = sqrt(693.0 / 128.0) * (-20.0 * x[3] * z[1] + 60.0 * x[1] * y[2] * z[1]);
                d2Ydydz[11] = sqrt(693.0 / 128.0) * (-20.0 * x[3] * y[1] + 20.0 * x[1] * y[3]);
                d2Ydzdz[11] = 0.0;
                d2Ydxdx[12] = sqrt(231.0 / 512.0) * (30.0 * x[4] - 180.0 * x[2] * y[2] + 30.0 * y[4]);
                d2Ydxdy[12] = sqrt(231.0 / 512.0) * (-120.0 * x[3] * y[1] + 120.0 * x[1] * y[3]);
                d2Ydxdz[12] = 0.0;
                d2Ydydy[12] = sqrt(231.0 / 512.0) * (-30.0 * x[4] + 180.0 * x[2] * y[2] - 30.0 * y[4]);
                d2Ydydz[12] = 0.0;
                d2Ydzdz[12] = 0.0;
              }
            }
            break;
          default:
            std::cout << "Angular momentum too high for transformation spherical harmonics." << std::endl;
            throw std::exception();
            break;
        }

        /* =======================
         *   Common Finalization
         * =======================*/
        for (unsigned int m = 0; m < 2 * angularMomentumOfMu + 1; ++m) {
          valueStore[m * offset] = radial * Y[m];
        }
        if (_highestDerivative >= 1) {
          for (unsigned int m = 0; m < 2 * angularMomentumOfMu + 1; ++m) {
            derivativeStore.x[m * offset] = radial * dYdx[m] + dradial * x[1] * Y[m];
            derivativeStore.y[m * offset] = radial * dYdy[m] + dradial * y[1] * Y[m];
            derivativeStore.z[m * offset] = radial * dYdz[m] + dradial * z[1] * Y[m];
          }
        }
        if (_highestDerivative >= 2) {
          for (unsigned int m = 0; m < 2 * angularMomentumOfMu + 1; ++m) {
            hessianStore.xx[m * offset] =
                radial * d2Ydxdx[m] + 2.0 * dradial * x[1] * dYdx[m] + ddradial * x[2] * Y[m] + dradial * Y[m];
            hessianStore.xy[m * offset] =
                radial * d2Ydxdy[m] + dradial * x[1] * dYdy[m] + dradial * y[1] * dYdx[m] + ddradial * x[1] * y[1] * Y[m];
            hessianStore.xz[m * offset] =
                radial * d2Ydxdz[m] + dradial * x[1] * dYdz[m] + dradial * z[1] * dYdx[m] + ddradial * x[1] * z[1] * Y[m];
            hessianStore.yy[m * offset] =
                radial * d2Ydydy[m] + 2.0 * dradial * y[1] * dYdy[m] + ddradial * y[2] * Y[m] + dradial * Y[m];
            hessianStore.yz[m * offset] =
                radial * d2Ydydz[m] + dradial * y[1] * dYdz[m] + dradial * z[1] * dYdy[m] + ddradial * y[1] * z[1] * Y[m];
            hessianStore.zz[m * offset] =
                radial * d2Ydzdz[m] + 2.0 * dradial * z[1] * dYdz[m] + ddradial * z[2] * Y[m] + dradial * Y[m];
          }
        }
      } /* if spherical*/
    }
  }
  thisBlockData->averageFunctionValues = thisBlockData->functionValues.cwiseAbs().colwise().mean();
  if (omp_get_thread_num() == 0)
    Timings::timeTaken("Tech. -    Basis On Grid Eval.");

  return thisBlockData;
}

} /* namespace Serenity */
