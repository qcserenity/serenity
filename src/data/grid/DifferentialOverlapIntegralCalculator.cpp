/**
 * @file DifferentialOverlapIntegralCalculator.cpp
 *
 * @date Jan 29, 2019
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
#include "data/grid/DifferentialOverlapIntegralCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"                   //   -- " --
#include "data/grid/BasisFunctionOnGridController.h" //Basis function values on a grid
#include "data/matrices/MatrixInBasis.h"             //ShellWiseAbsMax
#include "grid/GridController.h"                     //   -- " --
#include "misc/HelperFunctions.h"                    //Sparse map construction.

namespace Serenity {

void DifferentialOverlapIntegralCalculator::calculateDOI(const Eigen::MatrixXd& cX, const Eigen::MatrixXd& cY,
                                                         const Eigen::SparseMatrix<int>& basisFunctionToXMap,
                                                         const Eigen::SparseMatrix<int>& basisFunctionToYMap,
                                                         std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController,
                                                         Eigen::MatrixXd& dois) {
  takeTime("DOI -- calculations");
  const unsigned int nX = cX.cols();
  const unsigned int nY = cY.cols();
  assert(nX == basisFunctionToXMap.cols());
  assert(nY == basisFunctionToYMap.cols());

  const unsigned int nBlocks = basisFuncOnGridController->getNBlocks();
  const unsigned int nGridPts = basisFuncOnGridController->getNGridPoints();
  const auto& weights = basisFuncOnGridController->getGridController()->getWeights();
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> integralsPerThread(omp_get_max_threads(), Eigen::MatrixXd::Zero(nX, nY));
#pragma omp parallel for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    const unsigned int threadId = omp_get_thread_num();
    // Get block size.
    unsigned int blockEnd = 0;
    unsigned int firstOfBlock = basisFuncOnGridController->getFirstIndexOfBlock(block);
    if (block == nBlocks - 1) {
      blockEnd = nGridPts;
    }
    else {
      blockEnd = basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    const unsigned int nBlockPoints = blockEnd - firstOfBlock;
    const auto& basisFuncOnGrid = basisFuncOnGridController->getBlockOnGridData(block);
    /*
     * The idea behind this code is as follows:
     *   1. Extract the columns in x and y which have significant values on this block.
     *   2. Extract the dens-matrix containing the significant basis function values
     *      on this block.
     *   3. Contract significant x/y columns with  basis function values.
     *
     *   Note that this algorithm does not directly prescreen the coefficients of x and y
     *   on this block. The only information we use, is that there is at least one coefficient
     *   not removed. The basis function prescreening happens solely via the significant basis
     *   function values. Thus, the sparsity in the coefficients is not fully exploited. However,
     *   this allows us to write the entire contraction as matrix--matrix operations. Since this
     *   additional sparsity would only reduce the factor in front of the scaling law. I highly
     *   doubt that it is worth it to try to use it. I (MB) was not able to get the code faster
     *   than it is now.
     */

    // Significant basis function values on this block. Invert this presceening vector 1-->0 and 0-->1
    Eigen::VectorXi blockToBasisFunctionMap = basisFuncOnGrid->negligible;
    blockToBasisFunctionMap -= Eigen::VectorXi::Constant(blockToBasisFunctionMap.rows(), 1);
    blockToBasisFunctionMap = blockToBasisFunctionMap.array().abs();
    const Eigen::SparseVector<int> sparseBlockToBasisFunctionMap = blockToBasisFunctionMap.sparseView();

    // Significant x/y for this grid block.
    const Eigen::SparseVector<int> xToBlockMap = basisFunctionToXMap.transpose() * sparseBlockToBasisFunctionMap;
    const Eigen::SparseVector<int> yToBlockMap = basisFunctionToYMap.transpose() * sparseBlockToBasisFunctionMap;
    const Eigen::SparseMatrix<double> proXBlock = constructProjectionMatrixFromSparse(xToBlockMap);
    const Eigen::SparseMatrix<double> proYBlock = constructProjectionMatrixFromSparse(yToBlockMap);
    if (proXBlock.cols() < 1 || proYBlock.cols() < 1)
      continue;
    // Significant basis function values on this block.
    //    nFunc x nSigFunc
    const Eigen::SparseMatrix<double> proFunc = constructProjectionMatrixFromSparse(sparseBlockToBasisFunctionMap);
    const Eigen::MatrixXd basSig = basisFuncOnGrid->functionValues * proFunc; // nBlock x nSigFunc
    Eigen::MatrixXd xPoint = basSig * proFunc.transpose() * cX * proXBlock;   // nBlock x nSigX
    Eigen::MatrixXd yPoint = basSig * proFunc.transpose() * cY * proYBlock;   // nBlock x nSigY
    // Square
    xPoint.array() *= xPoint.array();
    yPoint.array() *= yPoint.array();
    // Contract with weights
    xPoint.array().colwise() *= weights.segment(firstOfBlock, nBlockPoints).array();
    // Calculate |X(r)|^2*|Y(r)|^2 * w(r) and add to matrix by back mapping to the original indices.
    integralsPerThread[threadId] += proXBlock * xPoint.transpose() * yPoint * proYBlock.transpose();
  } // for block
  dois = Eigen::MatrixXd::Zero(nX, nY);
  // sum over all threads
  for (unsigned int t = 0; t < integralsPerThread.size(); ++t) {
    dois += integralsPerThread[t];
  }
  dois = dois.array().sqrt();
  timeTaken(2, "DOI -- calculations");
}
void DifferentialOverlapIntegralCalculator::calculateDOI(std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController,
                                                         Eigen::MatrixXd& dois) {
  takeTime("DOI -- calculations");
  const unsigned int nBlocks = basisFuncOnGridController->getNBlocks();
  const unsigned int nGridPts = basisFuncOnGridController->getNGridPoints();
  const auto& weights = basisFuncOnGridController->getGridController()->getWeights();
  const unsigned int nX = basisFuncOnGridController->getNBasisFunctions();
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> intergralsPerThread(omp_get_max_threads(), Eigen::MatrixXd::Zero(nX, nX));
#pragma omp parallel for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    const unsigned int threadId = omp_get_thread_num();
    // Get block size.
    unsigned int blockEnd = 0;
    unsigned int firstOfBlock = basisFuncOnGridController->getFirstIndexOfBlock(block);
    if (block == nBlocks - 1) {
      blockEnd = nGridPts;
    }
    else {
      blockEnd = basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    const unsigned int nBlockPoints = blockEnd - firstOfBlock;
    const auto& basisFuncOnGrid = basisFuncOnGridController->getBlockOnGridData(block);
    /*
     * The idea behind this code is as follows:
     *   1. Extract the dens-matrix containing the significant basis function values
     *      on this block.
     *   2. Contract significant basis function values to get the DOIs for each basis function pair.
     */

    // Significant basis function values on this block. Invert this presceening vector 1-->0 and 0-->1
    Eigen::VectorXi blockToBasisFunctionMap = basisFuncOnGrid->negligible;
    blockToBasisFunctionMap -= Eigen::VectorXi::Constant(blockToBasisFunctionMap.rows(), 1);
    blockToBasisFunctionMap = blockToBasisFunctionMap.array().abs();
    const Eigen::SparseVector<int> sparseBlockToBasisFunctionMap = blockToBasisFunctionMap.sparseView();

    // Significant basis function values on this block.
    //    nFunc x nSigFunc
    const Eigen::SparseMatrix<double> proFunc = constructProjectionMatrixFromSparse(sparseBlockToBasisFunctionMap);
    Eigen::MatrixXd basSig = basisFuncOnGrid->functionValues * proFunc; // nBlock x nSigFunc
    // Square
    basSig.array() *= basSig.array();
    // Contract with weights
    Eigen::MatrixXd basSigWeighted = basSig;
    basSigWeighted.array().colwise() *= weights.segment(firstOfBlock, nBlockPoints).array();
    // Calculate |X(r)|^2*|Y(r)|^2 * w(r) and add to matrix by back mapping to the original indices.
    intergralsPerThread[threadId] += proFunc * basSigWeighted.transpose() * basSig * proFunc.transpose();
  } // for block
  dois = Eigen::MatrixXd::Zero(nX, nX);
  // sum over all threads
  for (unsigned int t = 0; t < intergralsPerThread.size(); ++t) {
    dois += intergralsPerThread[t];
  }
  dois = dois.array().sqrt();
  timeTaken(0, "DOI -- calculations");
}

std::vector<ShellPairData> DifferentialOverlapIntegralCalculator::calculateDOIShellPairData(
    std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController, double cutOff) {
  MatrixInBasis<RESTRICTED> dois(basisFuncOnGridController->getBasisController());
  calculateDOI(basisFuncOnGridController, dois);
  Eigen::MatrixXd shellWiseMax = dois.shellWiseAbsMax();
  const unsigned int nShells = basisFuncOnGridController->getBasisController()->getReducedNBasisFunctions();
  std::vector<ShellPairData> shellPairData;
  for (unsigned int iShell = 0; iShell < nShells; ++iShell) {
    for (unsigned int jShell = 0; jShell <= iShell; ++jShell) {
      if (shellWiseMax(iShell, jShell) > cutOff) {
        shellPairData.emplace_back(iShell, jShell, shellWiseMax(iShell, jShell));
      }
    }
  }
  return shellPairData;
}

} /* namespace Serenity */
