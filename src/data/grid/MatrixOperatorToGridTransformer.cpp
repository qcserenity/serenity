/**
 * @file   MatrixOperatorToGridTransformer.cpp
 * @author Thomas Dresselhaus, last rework Jan Unsleber
 *
 * @date   May 15, 2015 , last rework May 9, 2017
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
#include "data/grid/MatrixOperatorToGridTransformer.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "math/Matrix.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <cassert>
#include <cmath>

namespace Serenity {
/*********************************************************************************
 * Actual (private) calculation method                                           *
 * According to: O. Treutler and R. Ahlrichs, J. Chem. Phys (1995), 102, 346     *
 *********************************************************************************/
Eigen::SparseVector<int>
MatrixOperatorToGridTransformer::transform(const std::vector<std::reference_wrapper<const Eigen::MatrixXd>>& matrices,
                                           BasisFunctionOnGridController& basisFunctionOnGridController,
                                           std::vector<std::reference_wrapper<Eigen::VectorXd>>& results,
                                           std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>>* resultGradientsPtr,
                                           std::vector<Hessian<std::reference_wrapper<Eigen::VectorXd>>>* resultHessiansPtr) {
  /*
   * Some comments and variables below are written for the density, as this is what the referenced
   * paper is about. Don't be irritated by this.
   */
  const unsigned int nData = matrices.size();
  /*
   * Sanity checks
   */
  if (results.size() != nData)
    throw SerenityError("MatrixOperatorToGridTransformer: Mismatch of data field sizes.");
  if (resultGradientsPtr) {
    if (resultGradientsPtr->size() != nData)
      throw SerenityError("MatrixOperatorToGridTransformer: Mismatch of data field sizes.");
    if (basisFunctionOnGridController.getHighestDerivative() < 1 && !resultHessiansPtr)
      basisFunctionOnGridController.setHighestDerivative(1);
  }
  if (resultHessiansPtr) {
    if (resultHessiansPtr->size() != nData)
      throw SerenityError("MatrixOperatorToGridTransformer: Mismatch of data field sizes.");
    if (basisFunctionOnGridController.getHighestDerivative() < 2)
      basisFunctionOnGridController.setHighestDerivative(2);
  }

  /*
   * Clear/init output
   */
  const unsigned int nPoints = basisFunctionOnGridController.getNGridPoints();
  for (auto& result : results) {
    if (result.get().size() != nPoints)
      result.get().resize(nPoints);
    result.get().setZero();
  }
  if (resultGradientsPtr) {
    for (auto& grad : *resultGradientsPtr) {
      for (auto& component : grad) {
        if (component.get().size() != nPoints)
          component.get().resize(nPoints);
        component.get().setZero();
      }
    }
  }
  if (resultHessiansPtr) {
    for (auto& hess : *resultHessiansPtr) {
      for (auto& component : hess) {
        if (component.get().size() != nPoints)
          component.get().resize(nPoints);
        component.get().setZero();
      }
    }
  }
  /*
   * Loop over blocks of grid points
   */
  const unsigned int nBlocks = basisFunctionOnGridController.getNBlocks();
  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif
  std::vector<Eigen::Triplet<int>> newVec;
  std::vector<std::vector<Eigen::Triplet<int>>> nonNegligibleBlocks(nThreads, newVec);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int blockNumber = 0; blockNumber < nBlocks; ++blockNumber) {
    /*
     * Get in data for this block locally
     */
    const auto& thisBlockData = basisFunctionOnGridController.getBlockOnGridData(blockNumber);
    const unsigned int blockFirstIndex = basisFunctionOnGridController.getFirstIndexOfBlock(blockNumber);
    const unsigned int blockSize = thisBlockData->functionValues.rows();
    /*
     * This projection matrix projects to the non-negligible basis functions of any matrix expressed
     * in the associated basis. Note that this matrix is extremely sparse and multiplication with it
     * are extremely fast.
     */
    const Eigen::SparseMatrix<double> projection = constructProjectionMatrix(thisBlockData->negligible);
    if (projection.nonZeros() > 0) {
      unsigned int threadID = 0;
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      nonNegligibleBlocks[threadID].push_back(Eigen::Triplet<int>(blockNumber, 0, 1));
    }
    else {
      continue;
    }

    /*
     * Product of a basis function value on a grid point with the corresponding entries in the
     * density matrix.
     * Analog using the basis function gradients instead of the values. This is needed
     * for the (efficient) evaluation of the Hessian of the density.
     * (I.e. here is calculated: rho(point) * d/dx (nu(Point)) and analogues with d/dy and d/dz.)
     * See 'a' in the part 'A comment on the program structure' in:
     * O. Treutler and R. Ahlrichs, J. Chem. Phys (1995), 102, 346
     */
    const Eigen::MatrixXd sigFunctionValues = thisBlockData->functionValues * projection;
    Eigen::MatrixXd signGradX;
    Eigen::MatrixXd signGradY;
    Eigen::MatrixXd signGradZ;
    if (resultGradientsPtr || resultHessiansPtr) {
      signGradX = thisBlockData->derivativeValues->x * projection;
      signGradY = thisBlockData->derivativeValues->y * projection;
      signGradZ = thisBlockData->derivativeValues->z * projection;
    } // if resultGradientsPtr||resultHessiansPtr

    /*
     * Loop over data sets
     */
    for (unsigned int i = 0; i < nData; ++i) {
      /*
       * Abbreviate data to work on
       */
      const auto& thisMatrix = matrices[i].get();
      auto& thisResults = results[i].get();
      const Eigen::MatrixXd signMatrix = projection.transpose() * thisMatrix * projection;
      const Eigen::MatrixXd basis_P = sigFunctionValues * signMatrix;
      thisResults.segment(blockFirstIndex, blockSize) = (basis_P.array() * sigFunctionValues.array()).rowwise().sum();
      if (resultGradientsPtr) {
        auto& grad = (*resultGradientsPtr)[i];
        grad.x.get().segment(blockFirstIndex, blockSize) = 2.0 * (basis_P.array() * signGradX.array()).rowwise().sum();
        grad.y.get().segment(blockFirstIndex, blockSize) = 2.0 * (basis_P.array() * signGradY.array()).rowwise().sum();
        grad.z.get().segment(blockFirstIndex, blockSize) = 2.0 * (basis_P.array() * signGradZ.array()).rowwise().sum();
      } //  if resultGradientsPtr
      if (resultHessiansPtr) {
        auto& hess = (*resultHessiansPtr)[i];
        const Eigen::MatrixXd signHessXX = thisBlockData->secondDerivativeValues->xx * projection;
        const Eigen::MatrixXd signHessXY = thisBlockData->secondDerivativeValues->xy * projection;
        const Eigen::MatrixXd signHessXZ = thisBlockData->secondDerivativeValues->xz * projection;
        const Eigen::MatrixXd signHessYY = thisBlockData->secondDerivativeValues->yy * projection;
        const Eigen::MatrixXd signHessYZ = thisBlockData->secondDerivativeValues->yz * projection;
        const Eigen::MatrixXd signHessZZ = thisBlockData->secondDerivativeValues->zz * projection;
        const Eigen::MatrixXd grad_mat_x = signGradX * signMatrix;
        const Eigen::MatrixXd grad_mat_y = signGradY * signMatrix;
        const Eigen::MatrixXd grad_mat_z = signGradZ * signMatrix;
        hess.xx.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessXX.array() + grad_mat_x.array() * signGradX.array()).rowwise().sum();
        hess.xy.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessXY.array() + grad_mat_x.array() * signGradY.array()).rowwise().sum();
        hess.xz.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessXZ.array() + grad_mat_x.array() * signGradZ.array()).rowwise().sum();
        hess.yy.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessYY.array() + grad_mat_y.array() * signGradY.array()).rowwise().sum();
        hess.yz.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessYZ.array() + grad_mat_y.array() * signGradZ.array()).rowwise().sum();
        hess.zz.get().segment(blockFirstIndex, blockSize) =
            2.0 * (basis_P.array() * signHessZZ.array() + grad_mat_z.array() * signGradZ.array()).rowwise().sum();
      } // if resultHessiansPtr
    }   /*Data*/
  }     /*Block*/
  std::vector<Eigen::Triplet<int>> finalTripletSet = nonNegligibleBlocks[0];
  for (unsigned int iThread = 1; iThread < nonNegligibleBlocks.size(); ++iThread) {
    finalTripletSet.insert(finalTripletSet.end(), nonNegligibleBlocks[iThread].begin(), nonNegligibleBlocks[iThread].end());
  }
  Eigen::SparseMatrix<int> sparseBlocks(nBlocks, 1);
  sparseBlocks.setFromTriplets(finalTripletSet.begin(), finalTripletSet.end());
  return sparseBlocks.col(0);
}

} /* namespace Serenity */
