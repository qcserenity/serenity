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

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                                                                  const double blockAveThreshold)
  : _basisFunctionOnGridControllerA(basisFunctionOnGridController),
    _basisFunctionOnGridControllerB(basisFunctionOnGridController),
    _blockAveThreshold(blockAveThreshold) {
}

template<Options::SCF_MODES SCFMode>
ScalarOperatorToMatrixAdder<SCFMode>::ScalarOperatorToMatrixAdder(
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerA,
    std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerB, double blockAveThreshold)
  : _basisFunctionOnGridControllerA(basisFunctionOnGridControllerA),
    _basisFunctionOnGridControllerB(basisFunctionOnGridControllerB),
    _blockAveThreshold(blockAveThreshold) {
  if (!isDefinedOnSameGrid(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB))
    throw SerenityError("Both basis sets have to be evaluated on the same grid!");
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix,
                                                                     const GridPotential<SCFMode>& scalarOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  if (!isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA))
    throw SerenityError("The scalar operator and the BasisFunctionOnGridController passed to addScalarOperatorToMatrix "
                        "are not defined on the same grid!");

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  auto nbA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  auto nbB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  std::vector<SPMatrix<SCFMode>> threadMatrices(nThreads, SPMatrix<SCFMode>(nbA, nbB));

#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], scalarOperator);
  }
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
  if (!isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA))
    throw SerenityError("The scalar operator and the BasisFunctionOnGridController passed to addScalarOperatorToMatrix "
                        "are not defined on the same grid!");
  for (const auto& component : gradientOperator) {
    if (!isDefinedOnSameGrid(component, *_basisFunctionOnGridControllerA))
      throw SerenityError("The gradient operator and the BasisFunctionOnGridController passed to "
                          "addScalarOperatorToMatrix are not defined on the same grid!");
  }

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  auto nbA = _basisFunctionOnGridControllerA->getBasisController()->getNBasisFunctions();
  auto nbB = _basisFunctionOnGridControllerB->getBasisController()->getNBasisFunctions();
  std::vector<SPMatrix<SCFMode>> threadMatrices(nThreads, SPMatrix<SCFMode>(nbA, nbB));

#pragma omp parallel for schedule(dynamic)
  for (unsigned iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], scalarOperator, gradientOperator);
  }
  for (auto& mat : threadMatrices) {
    matrix += mat;
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToGradient(Eigen::MatrixXd& matrix,
                                                                       const Eigen::MatrixXd& atomBasisProjection,
                                                                       const DensityMatrix<SCFMode>& D,
                                                                       const GridPotential<SCFMode>& scalarOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  if (!isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA))
    throw SerenityError("The scalar operator and the BasisFunctionOnGridController passed to "
                        "addScalarOperatorToGradient are not defined on the same grid!");

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  std::vector<Eigen::MatrixXd> threadMatrices(nThreads, Eigen::MatrixXd::Zero(atomBasisProjection.rows(), 3));

#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addGradientBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], atomBasisProjection, D, scalarOperator);
  }
  for (auto& mat : threadMatrices) {
    matrix += mat;
  }
  Timings::timeTaken("Tech. -    Grid to Matrix Int.");
  return;
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addScalarOperatorToGradient(
    Eigen::MatrixXd& matrix, const Eigen::MatrixXd& atomBasisProjection, const DensityMatrix<SCFMode>& D,
    const GridPotential<SCFMode>& scalarOperator, const Gradient<GridPotential<SCFMode>>& gradientOperator) {
  Timings::takeTime("Tech. -    Grid to Matrix Int.");
  if (!isDefinedOnSameGrid(scalarOperator, *_basisFunctionOnGridControllerA))
    throw SerenityError("The scalar operator and the BasisFunctionOnGridController passed to "
                        "addScalarOperatorToGradient are not defined on the same grid!");
  for (const auto& component : gradientOperator) {
    if (!isDefinedOnSameGrid(component, *_basisFunctionOnGridControllerA))
      throw SerenityError("The gradient operator and the BasisFunctionOnGridController passed to "
                          "addScalarOperatorToGradient are not defined on the same grid!");
  }

  unsigned nBlocks = _basisFunctionOnGridControllerA->getNBlocks();
  unsigned nThreads = omp_get_max_threads();

  std::vector<Eigen::MatrixXd> threadMatrices(nThreads, Eigen::MatrixXd::Zero(atomBasisProjection.rows(), 3));

#pragma omp parallel for schedule(dynamic)
  for (unsigned iBlock = 0; iBlock < nBlocks; ++iBlock) {
    unsigned iThread = omp_get_thread_num();
    auto& blockDataA = _basisFunctionOnGridControllerA->getBlockOnGridData(iBlock);
    auto& blockDataB = _basisFunctionOnGridControllerB->getBlockOnGridData(iBlock);
    this->addGradientBlock(iBlock, blockDataA, blockDataB, threadMatrices[iThread], atomBasisProjection, D,
                           scalarOperator, gradientOperator);
  }
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
  // function values for each grid point/basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  // basis function negligibility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  // pA has dimensions nPoints * nNonNegligible
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

    // basisFuncVA has dimensions nPoints * nNonNegligible
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
  // function values for each grid point/basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  const auto& basisFunctionValuesB = blockDataB->functionValues;
  // gradient values for each grid point/basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& gradBasisFunctionValuesA = blockDataA->derivativeValues;
  const auto& gradBasisFunctionValuesB = blockDataB->derivativeValues;
  // basis function negligibility
  const auto& negligibleA = blockDataA->negligible;
  const auto& negligibleB = blockDataB->negligible;
  // pA has dimensions nPoints * nNonNegligible
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

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addGradientBlock(
    unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB, Eigen::MatrixXd& m_AB,
    const Eigen::MatrixXd& atomBasisProjection, const DensityMatrix<SCFMode>& D, const GridPotential<SCFMode>& scalarPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB);
  // Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  // function values for each grid point/basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  // @todo extend to properly handle two different bases
  (void)blockDataB;
  const auto& bfDerivativesA = *blockDataA->derivativeValues;
  // basis function negligibility
  const auto& negligibleA = blockDataA->negligible;
  // pA has dimensions nBasis * nNonNegligible
  const auto pA = constructProjectionMatrix(negligibleA);

  // number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  // Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  for_spin(scalarPart, D) {
    Eigen::VectorXd scalarW =
        scalarPart_spin.segment(iGridStart, blockSize).cwiseProduct(weights.segment(iGridStart, blockSize));
    double average = scalarW.cwiseAbs().sum() / blockSize;

    if (average < _blockAveThreshold)
      return;

    // basisFuncVA has dimensions nPoints * nNonNegligible
    Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
    // scalA has dimensions nNonNegligible * nPoints
    Eigen::MatrixXd scalA = basisFuncVA.transpose() * scalarW.asDiagonal();

    if (useSym) {
      Eigen::MatrixXd tmpx = bfDerivativesA.x * pA;
      Eigen::MatrixXd x = scalA * tmpx;
      Eigen::MatrixXd tmpy = bfDerivativesA.y * pA;
      Eigen::MatrixXd y = scalA * tmpy;
      Eigen::MatrixXd tmpz = bfDerivativesA.z * pA;
      Eigen::MatrixXd z = scalA * tmpz;
      m_AB.col(0).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * x * pA.transpose()).colwise().sum()).transpose();
      m_AB.col(1).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * y * pA.transpose()).colwise().sum()).transpose();
      m_AB.col(2).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * z * pA.transpose()).colwise().sum()).transpose();
    }
    // else {
    // case of two different bases goes here
    // }
  }; /* for_spin */
}

template<Options::SCF_MODES SCFMode>
void ScalarOperatorToMatrixAdder<SCFMode>::addGradientBlock(
    unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
    std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB, Eigen::MatrixXd& m_AB,
    const Eigen::MatrixXd& atomBasisProjection, const DensityMatrix<SCFMode>& D,
    const GridPotential<SCFMode>& scalarPart, const Gradient<GridPotential<SCFMode>>& gradientPart) {
  bool useSym = isDefinedInSameBasis(*_basisFunctionOnGridControllerA, *_basisFunctionOnGridControllerB);
  // Get weights
  const auto& weights = _basisFunctionOnGridControllerA->getGridController()->getWeights();
  // function values for each grid point/basis function combination (Dimension: nPoints x nBasisFunctions)
  const auto& basisFunctionValuesA = blockDataA->functionValues;
  // @todo extend to properly handle two different bases
  (void)blockDataB;
  // gradient values for each grid point/basis function combination (Dimension: 3 * nPoints x nBasisFunctions)
  const auto& bfDerivativesA = *blockDataA->derivativeValues;
  const Hessian<Eigen::MatrixXd>& bfSecondDer = *blockDataA->secondDerivativeValues;
  // basis function negligibility
  const auto& negligibleA = blockDataA->negligible;
  // pA has dimensions nPoints * nNonNegligible
  const auto pA = constructProjectionMatrix(negligibleA);
  // number of grid points in this block
  const unsigned int blockSize = blockDataA->functionValues.rows();
  // Get first index of this Block
  const unsigned int iGridStart = _basisFunctionOnGridControllerA->getFirstIndexOfBlock(iBlock);
  const auto& gradientPartX = gradientPart.x;
  const auto& gradientPartY = gradientPart.y;
  const auto& gradientPartZ = gradientPart.z;
  for_spin(scalarPart, gradientPartX, gradientPartY, gradientPartZ, D) {
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

    Eigen::MatrixXd basisFuncVA = basisFunctionValuesA * pA;
    Eigen::MatrixXd gradBasisFuncVA_x = bfDerivativesA.x * pA;
    Eigen::MatrixXd gradBasisFuncVA_y = bfDerivativesA.y * pA;
    Eigen::MatrixXd gradBasisFuncVA_z = bfDerivativesA.z * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_xx = bfSecondDer.xx * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_xy = bfSecondDer.xy * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_xz = bfSecondDer.xz * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_yy = bfSecondDer.yy * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_yz = bfSecondDer.yz * pA;
    Eigen::MatrixXd secondGradBasisFuncVA_zz = bfSecondDer.zz * pA;

    Eigen::MatrixXd scalA =
        basisFuncVA.transpose() * scalarW.asDiagonal() + gradBasisFuncVA_x.transpose() * gradWx.asDiagonal() +
        gradBasisFuncVA_y.transpose() * gradWy.asDiagonal() + gradBasisFuncVA_z.transpose() * gradWz.asDiagonal();
    Eigen::MatrixXd gradA_x = basisFuncVA.transpose() * gradWx.asDiagonal();
    Eigen::MatrixXd gradA_y = basisFuncVA.transpose() * gradWy.asDiagonal();
    Eigen::MatrixXd gradA_z = basisFuncVA.transpose() * gradWz.asDiagonal();

    if (useSym) {
      Eigen::MatrixXd x = scalA * gradBasisFuncVA_x + gradA_x * secondGradBasisFuncVA_xx +
                          gradA_y * secondGradBasisFuncVA_xy + gradA_z * secondGradBasisFuncVA_xz;
      Eigen::MatrixXd y = scalA * gradBasisFuncVA_y + gradA_x * secondGradBasisFuncVA_xy +
                          gradA_y * secondGradBasisFuncVA_yy + gradA_z * secondGradBasisFuncVA_yz;
      Eigen::MatrixXd z = scalA * gradBasisFuncVA_z + gradA_x * secondGradBasisFuncVA_xz +
                          gradA_y * secondGradBasisFuncVA_yz + gradA_z * secondGradBasisFuncVA_zz;

      m_AB.col(0).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * x * pA.transpose()).colwise().sum()).transpose();
      m_AB.col(1).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * y * pA.transpose()).colwise().sum()).transpose();
      m_AB.col(2).noalias() -=
          atomBasisProjection * (D_spin.cwiseProduct(pA * z * pA.transpose()).colwise().sum()).transpose();
    }
    // else {
    // case of two different bases goes here
    // }
    // }
  }; /* for_spin */
}

template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::RESTRICTED>;
template class ScalarOperatorToMatrixAdder<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
