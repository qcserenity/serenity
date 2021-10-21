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
#ifdef _OPENMP
  omp_init_lock(&_lock);
#endif
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
#ifdef _OPENMP
  omp_init_lock(&_lock);
#endif
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix,
                                                                     const GridPotential<SCFMode>& scalarOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  // checks
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));
  // Get number of grid blocks
  const unsigned int nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  // Build matrix for parallel computation
  unsigned int nThreads = (unsigned int)omp_get_max_threads();
  std::vector<std::unique_ptr<SPMatrix<SCFMode>>> threadMatrices(nThreads);

#ifdef _OPENMP
  Eigen::setNbThreads(1);
#endif
#pragma omp parallel
  {
    unsigned int threadID = omp_get_thread_num();
    threadMatrices[threadID].reset(
        new SPMatrix<SCFMode>(_basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions(),
                              _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions()));

    // Loop over blocks of grid points
#pragma omp for schedule(static, 1)
    for (unsigned int blockNumber = 0; blockNumber < nBlocks; ++blockNumber) {
      // Data for this block
      auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(blockNumber);
      auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(blockNumber);
      // Integrate and add to Fock matrix for this thread
      addBlock(blockNumber, blockDataA, blockDataB, *threadMatrices[threadID], scalarOperator);
    } /* loop over all blocks of grid points */
  }   // end pragma omp parallel
#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif
  for (auto& mat : threadMatrices) {
    if (mat != nullptr) {
      matrix += *mat;
    }
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(
    // TODO maybe a different type?
    SPMatrix<SCFMode>& matrix, const GridPotential<SCFMode>& scalarOperator,
    const Gradient<GridPotential<SCFMode>>& gradientOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  // checks
  assert(scalarOperator.isValid());
  assert(isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA));
  assert(_basisFunctionOnGridControllerA->getHighestDerivative() >= 1);
  assert(_basisFunctionOnGridControllerB->getHighestDerivative() >= 1);
  for (const auto& component : gradientOperator) {
    (void)component;
    assert(component.isValid());
    assert(isDefinedOnSameGrid(component, *_basisFunctionOnGridControllerA));
  }
  // Get number of grid blocks
  const unsigned int nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  // Build matrix for parallel computation
  unsigned int nThreads = (unsigned int)omp_get_max_threads();
  std::vector<std::unique_ptr<SPMatrix<SCFMode>>> threadMatrices(nThreads);

#ifdef _OPENMP
  Eigen::setNbThreads(1);
#endif
#pragma omp parallel
  {
    unsigned int threadID = omp_get_thread_num();
    threadMatrices[threadID].reset(
        new SPMatrix<SCFMode>(_basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions(),
                              _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions()));

    // Loop over blocks of grid points
#pragma omp for schedule(static, 1)
    for (unsigned int blockNumber = 0; blockNumber < nBlocks; ++blockNumber) {
      // Data for this block
      auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(blockNumber);
      auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(blockNumber);
      // Integrate and add to Fock matrix for this thread
      addBlock(blockNumber, blockDataA, blockDataB, *threadMatrices[threadID], scalarOperator, gradientOperator);
    } /* loop over all blocks of grid points */
  }   // end pragma omp parallel
#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif
  for (auto& mat : threadMatrices) {
    if (mat != nullptr) {
      matrix += *mat;
    }
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
    // Pre-multiply operator and weights
    Eigen::VectorXd scalarW(blockSize);
    scalarW.array() = scalarPart_spin.segment(iGridStart, blockSize).array() * weights.segment(iGridStart, blockSize).array();
    double average = scalarW.cwiseAbs().sum();
    average /= blockSize;
    // cycle for_spin macro if non-significant
    if (average < _blockAveThreshold)
      return;
    // Contract with AO basis functions
    if (!useSym) {
      const Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
      const Eigen::MatrixXd basisFuncVB = basisFunctionValuesB * pB;
      const Eigen::MatrixXd tmp = (basisFuncVA.array().colwise() * scalarW.array()).matrix().transpose() * basisFuncVB;
      m_AB_spin += pA * tmp * pB.transpose();
    }
    else {
      const Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
      const Eigen::MatrixXd tmp = (basisFuncVA.array().colwise() * scalarW.array()).matrix().transpose() * basisFuncVA;
      m_AB_spin += pA * tmp * pB.transpose();
    }
  }; /* for_spin */
  return;
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
    Eigen::VectorXd scalarW(blockSize);
    auto gradW = makeGradient<Eigen::VectorXd>(blockSize);
    scalarW.array() = scalarPart_spin.segment(iGridStart, blockSize).array() * weights.segment(iGridStart, blockSize).array();
    gradW.x.array() =
        gradientPartX_spin.segment(iGridStart, blockSize).array() * weights.segment(iGridStart, blockSize).array();
    gradW.y.array() =
        gradientPartY_spin.segment(iGridStart, blockSize).array() * weights.segment(iGridStart, blockSize).array();
    gradW.z.array() =
        gradientPartZ_spin.segment(iGridStart, blockSize).array() * weights.segment(iGridStart, blockSize).array();
    // Calculate average for prescreening.
    double average = scalarW.cwiseAbs().sum();
    average += gradW.x.cwiseAbs().sum();
    average += gradW.y.cwiseAbs().sum();
    average += gradW.z.cwiseAbs().sum();
    average /= blockSize;
    // cycle for_spin macro if non-significant
    if (average < _blockAveThreshold)
      return;

    Eigen::MatrixXd grad_AB;
    Eigen::MatrixXd scalar_AB;
    if (useSym) {
      const Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
      const Eigen::MatrixXd gradBasisFuncVA_x = gradBasisFunctionValuesA->x * pA;
      const Eigen::MatrixXd gradBasisFuncVA_y = gradBasisFunctionValuesA->y * pA;
      const Eigen::MatrixXd gradBasisFuncVA_z = gradBasisFunctionValuesA->z * pA;
      scalar_AB = (basisFuncVA.array().colwise() * scalarW.array()).matrix().transpose() * basisFuncVA;
      const Eigen::MatrixXd grad_A = gradBasisFuncVA_x.array().colwise() * gradW.x.array() +
                                     gradBasisFuncVA_y.array().colwise() * gradW.y.array() +
                                     gradBasisFuncVA_z.array().colwise() * gradW.z.array();
      grad_AB = basisFuncVA.transpose() * grad_A + grad_A.transpose() * basisFuncVA;
    }
    else {
      const Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
      const Eigen::MatrixXd basisFuncVB = basisFunctionValuesB * pB;
      const Eigen::MatrixXd gradBasisFuncVA_x = gradBasisFunctionValuesA->x * pA;
      const Eigen::MatrixXd gradBasisFuncVA_y = gradBasisFunctionValuesA->y * pA;
      const Eigen::MatrixXd gradBasisFuncVA_z = gradBasisFunctionValuesA->z * pA;
      const Eigen::MatrixXd gradBasisFuncVB_x = gradBasisFunctionValuesB->x * pB;
      const Eigen::MatrixXd gradBasisFuncVB_y = gradBasisFunctionValuesB->y * pB;
      const Eigen::MatrixXd gradBasisFuncVB_z = gradBasisFunctionValuesB->z * pB;
      scalar_AB = (basisFuncVA.array().colwise() * scalarW.array()).matrix().transpose() * basisFuncVB;
      const Eigen::MatrixXd grad_A = gradBasisFuncVA_x.array().colwise() * gradW.x.array() +
                                     gradBasisFuncVA_y.array().colwise() * gradW.y.array() +
                                     gradBasisFuncVA_z.array().colwise() * gradW.z.array();
      const Eigen::MatrixXd grad_B = gradBasisFuncVB_x.array().colwise() * gradW.x.array() +
                                     gradBasisFuncVB_y.array().colwise() * gradW.y.array() +
                                     gradBasisFuncVB_z.array().colwise() * gradW.z.array();
      grad_AB = basisFuncVA.transpose() * grad_B + grad_A.transpose() * basisFuncVB;
    }
    m_AB_spin += pA * (scalar_AB + grad_AB) * pB.transpose();
  }; /* for_spin */
  return;
}

template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED>;
template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
