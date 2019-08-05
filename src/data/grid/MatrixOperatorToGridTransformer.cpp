/**
 * @file   MatrixOperatorToGridTransformer.cpp
 * @author Thomas Dresselhaus, last rework Jan Unsleber
 *
 * @date   May 15, 2015 , last rework May 9, 2017
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
#include "data/grid/MatrixOperatorToGridTransformer.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "math/Matrix.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <cassert>
#include <cmath>

namespace Serenity {
/*********************************************************************************
 * Actual (private) calculation method                                           *
 * According to: O. Treutler and R. Ahlrichs, J. Chem. Phys (1995), 102, 346     *
 *********************************************************************************/
void MatrixOperatorToGridTransformer::transform(const std::vector<std::reference_wrapper<const Eigen::MatrixXd>>& matrices,
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
  }
  if (resultGradientsPtr) {
    for (auto& grad : *resultGradientsPtr) {
      for (auto& component : grad) {
        if (component.get().size() != nPoints)
          component.get().resize(nPoints);
      }
    }
  }
  if (resultHessiansPtr) {
    for (auto& hess : *resultHessiansPtr) {
      for (auto& component : hess) {
        if (component.get().size() != nPoints)
          component.get().resize(nPoints);
      }
    }
  }
  const unsigned int nBasisFunctions = basisFunctionOnGridController.getNBasisFunctions();
  /*
   * Loop over blocks of grid points
   */
  const unsigned int nBlocks = basisFunctionOnGridController.getNBlocks();
  ;
#pragma omp parallel for schedule(dynamic)
  for (unsigned int blockNumber = 0; blockNumber < nBlocks; ++blockNumber) {
    /*
     * Get in data for this block locally
     */
    const auto& thisBlockData = basisFunctionOnGridController.getBlockOnGridData(blockNumber);
    const unsigned int blockFirstIndex = basisFunctionOnGridController.getFirstIndexOfBlock(blockNumber);
    const unsigned int blockSize = thisBlockData->functionValues.rows();

    assert(thisBlockData->functionValues.cols() == nBasisFunctions);
    /*
     * Product of a basis function value on a grid point with the corresponding entries in the
     * density matrix.
     * Analog using the basis function gradients instead of the values. This is needed
     * for the (efficient) evaluation of the Hessian of the density.
     * (I.e. here is calculated: rho(point) * d/dx (nu(Point)) and analogues with d/dy and d/dz.)
     * See 'a' in the part 'A comment on the program structure' in:
     * O. Treutler and R. Ahlrichs, J. Chem. Phys (1995), 102, 346
     */
    Eigen::MatrixXd nuAndP(blockSize, (resultHessiansPtr ? 4 : 1) * nBasisFunctions);

    /*
     * Loop over data sets
     */
    for (unsigned int i = 0; i < nData; ++i) {
      /*
       * Abbreviate data to work on
       */
      const auto& thisMatrix = matrices[i].get();
      auto& thisResults = results[i].get();

      // (Re-)init nuAndP
      nuAndP.setZero();

      /*====================================
       *  Pre-calculations incl. screening
       *====================================*/
      for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
        // screening
        if (thisBlockData->negligible[nu])
          continue;
        for (unsigned int mu = 0; mu < nBasisFunctions; ++mu) {
          const double rhoMuNu = thisMatrix.data()[mu * nBasisFunctions + nu];
          nuAndP.col(mu) += rhoMuNu * thisBlockData->functionValues.col(nu);

          if (resultHessiansPtr) {
            nuAndP.col(mu + 1 * nBasisFunctions) += rhoMuNu * thisBlockData->derivativeValues->x.col(nu);
            nuAndP.col(mu + 2 * nBasisFunctions) += rhoMuNu * thisBlockData->derivativeValues->y.col(nu);
            nuAndP.col(mu + 3 * nBasisFunctions) += rhoMuNu * thisBlockData->derivativeValues->z.col(nu);
          }
        }
      }

      /*==========
       *  Values
       *==========*/
      thisResults.segment(blockFirstIndex, blockSize).setZero();
      for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
        if (thisBlockData->negligible[nu])
          continue;
        thisResults.segment(blockFirstIndex, blockSize).array() +=
            thisBlockData->functionValues.col(nu).array() * nuAndP.col(nu).array();
      }

      /*====================
       *  First Derivatives
       *====================*/
      if (resultGradientsPtr) {
        auto& grad = (*resultGradientsPtr)[i];
        grad.x.get().segment(blockFirstIndex, blockSize).setZero();
        grad.y.get().segment(blockFirstIndex, blockSize).setZero();
        grad.z.get().segment(blockFirstIndex, blockSize).setZero();
        for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
          if (thisBlockData->negligible[nu])
            continue;
          grad.x.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * thisBlockData->derivativeValues->x.col(nu).array() * nuAndP.col(nu).array();
          grad.y.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * thisBlockData->derivativeValues->y.col(nu).array() * nuAndP.col(nu).array();
          grad.z.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * thisBlockData->derivativeValues->z.col(nu).array() * nuAndP.col(nu).array();
        }
      }
      /*=====================
       *  Second Derivative
       *=====================*/
      if (resultHessiansPtr) {
        auto& hess = (*resultHessiansPtr)[i];
        hess.xx.get().segment(blockFirstIndex, blockSize).setZero();
        hess.xy.get().segment(blockFirstIndex, blockSize).setZero();
        hess.xz.get().segment(blockFirstIndex, blockSize).setZero();
        hess.yy.get().segment(blockFirstIndex, blockSize).setZero();
        hess.yz.get().segment(blockFirstIndex, blockSize).setZero();
        hess.zz.get().segment(blockFirstIndex, blockSize).setZero();
        for (unsigned int nu = 0; nu < nBasisFunctions; ++nu) {
          if (thisBlockData->negligible[nu])
            continue;
          hess.xx.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->xx.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->x.col(nu).array() * nuAndP.col(1 * nBasisFunctions + nu).array());
          hess.xy.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->xy.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->y.col(nu).array() * nuAndP.col(1 * nBasisFunctions + nu).array());
          hess.xz.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->xz.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->z.col(nu).array() * nuAndP.col(1 * nBasisFunctions + nu).array());
          hess.yy.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->yy.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->y.col(nu).array() * nuAndP.col(2 * nBasisFunctions + nu).array());
          hess.yz.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->yz.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->z.col(nu).array() * nuAndP.col(2 * nBasisFunctions + nu).array());
          hess.zz.get().segment(blockFirstIndex, blockSize).array() +=
              2.0 * (thisBlockData->secondDerivativeValues->zz.col(nu).array() * nuAndP.col(nu).array() +
                     thisBlockData->derivativeValues->z.col(nu).array() * nuAndP.col(3 * nBasisFunctions + nu).array());
        }
      }
    } /*Data*/
  }   /*Block*/
}

} /* namespace Serenity */
