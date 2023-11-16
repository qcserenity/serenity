/**
 * @file   ScalarOperatorToMatrixAdder.cpp
 *
 * @date   Mar 22, 2014, Apr 20, 2017
 * @author Thomas Dresselhaus, M. Boeckers
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
#include "data/grid/ScalarOperatorToMatrixAdder.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/matrices/FockMatrix.h"
#include "grid/GridController.h"
#include "math/FloatMaths.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                                  const double blockAveThreshold)
  : _basisFunctionOnGridControllerA(basisFunctionOnGridController),
    _basisFunctionOnGridControllerB(basisFunctionOnGridController),
    _blockAveThreshold(blockAveThreshold) {
  assert(_basisFunctionOnGridControllerA);
}

template<Options::SCF_MODES SCFMode>
ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerA,
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerB, double blockAveThreshold)
  : _basisFunctionOnGridControllerA(basisFunctionOnGridControllerA),
    _basisFunctionOnGridControllerB(basisFunctionOnGridControllerB),
    _blockAveThreshold(blockAveThreshold) {
  assert(_basisFunctionOnGridControllerA);
  assert(_basisFunctionOnGridControllerB);
  assert(isDefinedOnSameGrid(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB) &&
         "Both basis sets have to be evaluated on the same grid!");
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix,
                                                                     const GridPotential<SCFMode>& scalarOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  auto nbA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  auto nbB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  std::vector<SPMatrix<SCFMode>> threadMatrices(nThreads, SPMatrix<SCFMode>(nbA, nbB));

  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], scalarOperator);
  }
  Eigen::setNbThreads(0);
  for (auto& mat : threadMatrices) {
    matrix += mat;
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix,
                                                                     const GridPotential<SCFMode>& scalarOperator,
                                                                     const Gradient<GridPotential<SCFMode>>& gradientOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));
  assert(_basisFunctionOnGridControllerA->getHighestDerivative() >= 1);
  assert(_basisFunctionOnGridControllerB->getHighestDerivative() >= 1);
  for (const auto& component : gradientOperator) {
    (void)component;
    assert(component.isValid());
    assert(isDefinedOnSameGrid(component, *_basisFunctionOnGridControllerA));
  }

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  auto nbA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  auto nbB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  std::vector<SPMatrix<SCFMode>> threadMatrices(nThreads, SPMatrix<SCFMode>(nbA, nbB));

  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], scalarOperator, gradientOperator);
  }
  Eigen::setNbThreads(0);
  for (auto& mat : threadMatrices) {
    matrix += mat;
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
    unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB, SPMatrix<SCFMode>& m_AB,
    const GridPotential<SCFMode>& scalarPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB);
  // Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  // function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  // basis function negligebility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  const auto pA = constructProjectionMatrix(negligibleA);
  const auto pB = constructProjectionMatrix(negligibleB);

  // number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  // Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  for_spin(m_AB, scalarPart) {
    Eigen::VectorXd scalarW =
        scalarPart_spin.segment(iGridStart, blockSize).cwiseProduct(weights.segment(iGridStart, blockSize));
    double average = scalarW.cwiseAbs().sum() / blockSize;

    if (average < _blockAveThreshold)
      return;

    Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
    Eigen::MatrixXd scalA = basisFuncVA.transpose() * scalarW.asDiagonal();

    if (useSym) {
      Eigen::MatrixXd tmp = scalA * basisFuncVA;
      m_AB_spin.noalias() += pA * tmp * pB.transpose();
    }
    else {
      Eigen::MatrixXd basisFuncVB = basisFunctionValuesB * pB;
      Eigen::MatrixXd tmp = scalA * basisFuncVB;
      m_AB_spin.noalias() += pA * tmp * pB.transpose();
    }
  }; /* for_spin */
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addBlock(
    unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB, SPMatrix<SCFMode>& m_AB,
    const GridPotential<SCFMode>& scalarPart, const Gradient<GridPotential<SCFMode>>& gradientPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB);
  // Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  // function values for each grid point/ basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  // gradient values for each grid point/ basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValuesA = blockDataA->derivativeValues;
  const auto& gradBasisFunctionValuesB = blockDataB->derivativeValues;
  // basis function negligebility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  const auto pA = constructProjectionMatrix(negligibleA);
  const auto pB = constructProjectionMatrix(negligibleB);
  // number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  // Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  auto& gradientPartX = gradientPart.x;
  auto& gradientPartY = gradientPart.y;
  auto& gradientPartZ = gradientPart.z;
  for_spin(m_AB, scalarPart, gradientPartX, gradientPartY, gradientPartZ) {
    // Pre-multiply operator and weights
    Eigen::VectorXd scalarW =
        weights.segment(iGridStart, blockSize).cwiseProduct(scalarPart_spin.segment(iGridStart, blockSize));
    Eigen::VectorXd gradWx =
        weights.segment(iGridStart, blockSize).cwiseProduct(gradientPartX_spin.segment(iGridStart, blockSize));
    Eigen::VectorXd gradWy =
        weights.segment(iGridStart, blockSize).cwiseProduct(gradientPartY_spin.segment(iGridStart, blockSize));
    Eigen::VectorXd gradWz =
        weights.segment(iGridStart, blockSize).cwiseProduct(gradientPartZ_spin.segment(iGridStart, blockSize));
    // Calculate average for prescreening.
    double average = scalarW.cwiseAbs().sum();
    average += gradWx.cwiseAbs().sum();
    average += gradWy.cwiseAbs().sum();
    average += gradWz.cwiseAbs().sum();
    // cycle for_spin macro if non-significant
    if (average / blockSize < _blockAveThreshold)
      return;

    Eigen::MatrixXd scalar_grad_AB;

    Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
    Eigen::MatrixXd gradBasisFuncVA_x = gradBasisFunctionValuesA->x * pA;
    Eigen::MatrixXd gradBasisFuncVA_y = gradBasisFunctionValuesA->y * pA;
    Eigen::MatrixXd gradBasisFuncVA_z = gradBasisFunctionValuesA->z * pA;
    Eigen::MatrixXd grad_A = gradWx.asDiagonal() * gradBasisFuncVA_x;
    grad_A.noalias() += gradWy.asDiagonal() * gradBasisFuncVA_y;
    grad_A.noalias() += gradWz.asDiagonal() * gradBasisFuncVA_z;

    if (useSym) {
      grad_A.noalias() += 0.5 * scalarW.asDiagonal() * basisFuncVA;
      Eigen::MatrixXd tmp = basisFuncVA.transpose() * grad_A;
      scalar_grad_AB = tmp;
      scalar_grad_AB += tmp.transpose();
    }
    else {
      Eigen::MatrixXd basisFuncVB = basisFunctionValuesB * pB;
      Eigen::MatrixXd gradBasisFuncVB_x = gradBasisFunctionValuesB->x * pB;
      Eigen::MatrixXd gradBasisFuncVB_y = gradBasisFunctionValuesB->y * pB;
      Eigen::MatrixXd gradBasisFuncVB_z = gradBasisFunctionValuesB->z * pB;

      scalar_grad_AB = basisFuncVA.transpose() * scalarW.asDiagonal() * basisFuncVB;

      Eigen::MatrixXd grad_B = gradWx.asDiagonal() * gradBasisFuncVB_x;
      grad_B.noalias() += gradWy.asDiagonal() * gradBasisFuncVB_y;
      grad_B.noalias() += gradWz.asDiagonal() * gradBasisFuncVB_z;

      scalar_grad_AB.noalias() += basisFuncVA.transpose() * grad_B;
      scalar_grad_AB.noalias() += grad_A.transpose() * basisFuncVB;
    }
    m_AB_spin.noalias() += pA * scalar_grad_AB * pB.transpose();
  }; /* for_spin */
}

template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED>;
template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
