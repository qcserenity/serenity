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
#include "grid/GridController.h"                     //   -- " --

namespace Serenity {

void DifferentialOverlapIntegralCalculator::calculateDOI(const Eigen::MatrixXd& c_x, const Eigen::MatrixXd& c_y,
                                                         const Eigen::SparseMatrix<int>& basisFunctionToXMap,
                                                         const Eigen::SparseMatrix<int>& basisFunctionToYMap,
                                                         std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController,
                                                         Eigen::MatrixXd& dois) {
  takeTime("DOI -- calculations");
  const unsigned int nX = c_x.cols();
  const unsigned int nY = c_y.cols();
  assert(nX == basisFunctionToXMap.cols());
  assert(nY == basisFunctionToYMap.cols());

  const unsigned int nBlocks = basisFuncOnGridController->getNBlocks();
  const unsigned int nGridPts = basisFuncOnGridController->getNGridPoints();
  const auto& weights = basisFuncOnGridController->getGridController()->getWeights();
  const unsigned int nBasisFunctions = basisFuncOnGridController->getBasisController()->getNBasisFunctions();

#ifdef _OPENMP
  // create a vector of matrices for each thread
  std::vector<Eigen::MatrixXd> z_lms(omp_get_max_threads(), Eigen::MatrixXd::Zero(nX, nY));
  unsigned int nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#else
  // or just one
  std::vector<Eigen::MatrixXd> z_lms(1, Eigen::MatrixXd::Zero(nX, nY));
#endif
#pragma omp parallel for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    const auto& basisFuncOnGrid = basisFuncOnGridController->getBlockOnGridData(block);
    // block --> phi
    Eigen::VectorXi blockToBasisFunctionMap = basisFuncOnGrid->negligible;
    // tmp
    blockToBasisFunctionMap -= Eigen::VectorXi::Constant(blockToBasisFunctionMap.rows(), 1);
    blockToBasisFunctionMap = blockToBasisFunctionMap.array().abs();
    // Calculate the intersection of the block-->Basis Function and x/y-->Basis Function maps.
    Eigen::SparseMatrix<int> phiToXMap_block(nBasisFunctions, nX);
    Eigen::SparseMatrix<int> phiToYMap_block(nBasisFunctions, nY);

    for (unsigned int x = 0; x < c_x.cols(); ++x) {
      phiToXMap_block.col(x) = Eigen::VectorXi(basisFunctionToXMap.col(x)).cwiseProduct(blockToBasisFunctionMap).sparseView();
    }
    for (unsigned int y = 0; y < c_y.cols(); ++y) {
      phiToYMap_block.col(y) = Eigen::VectorXi(basisFunctionToYMap.col(y)).cwiseProduct(blockToBasisFunctionMap).sparseView();
    }
    Eigen::SparseVector<int> xToBlockMap = basisFunctionToXMap.transpose().eval() * blockToBasisFunctionMap.sparseView();
    Eigen::SparseVector<int> yToBlockMap = basisFunctionToYMap.transpose().eval() * blockToBasisFunctionMap.sparseView();
    // Loop over block
    unsigned int blockEnd = 0;
    unsigned int firstOfBlock = basisFuncOnGridController->getFirstIndexOfBlock(block);
    if (block == nBlocks - 1) {
      blockEnd = nGridPts;
    }
    else {
      blockEnd = basisFuncOnGridController->getFirstIndexOfBlock(block + 1);
    }
    const unsigned int nBlockPoints = blockEnd - firstOfBlock;

    Eigen::MatrixXd x_point = Eigen::MatrixXd::Zero(nBlockPoints, nX);
    Eigen::MatrixXd y_point = Eigen::MatrixXd::Zero(nBlockPoints, nY);
    // Loop over x in map block-->x
    for (Eigen::SparseVector<int>::InnerIterator itX(xToBlockMap); itX; ++itX) {
      unsigned int x = itX.row();
      // Loop over phi in map x--> phi intersection block-->phi
      for (Eigen::SparseMatrix<int>::InnerIterator itPhi(phiToXMap_block, x); itPhi; ++itPhi) {
        x_point.col(x) += c_x(itPhi.row(), x) * basisFuncOnGrid->functionValues.col(itPhi.row()).eval();
      } // for itPhi
    }   // for itX
    // Loop over y in map block-->Y
    for (Eigen::SparseVector<int>::InnerIterator itY(yToBlockMap); itY; ++itY) {
      unsigned int y = itY.row();
      // Loop over phi in map y--> phi intersection block-->phi
      for (Eigen::SparseMatrix<int>::InnerIterator itPhi(phiToYMap_block, y); itPhi; ++itPhi) {
        y_point.col(y) += c_y(itPhi.row(), y) * basisFuncOnGrid->functionValues.col(itPhi.row()).eval();
      } // for itX
    }
    const Eigen::VectorXd blockWeights = weights.segment(firstOfBlock, nBlockPoints).eval();
    // Add |X(r)|^2*|Y(r)|^2 * w(r) to square matrix.
    // Loop over all x in map block--> x and all y in map block-->y
    for (Eigen::SparseVector<int>::InnerIterator itX(xToBlockMap); itX; ++itX) {
      // x(r)*w(r)*x(r)
      unsigned int x = itX.row();
      const Eigen::VectorXd xSquare = x_point.col(x).cwiseProduct(blockWeights).cwiseProduct(x_point.col(x)).eval();
      for (Eigen::SparseVector<int>::InnerIterator itY(yToBlockMap); itY; ++itY) {
        // y(r)*y(r)
        unsigned int y = itY.row();
        const Eigen::VectorXd ySquare = y_point.col(y).cwiseProduct(y_point.col(y)).eval();
        // z(x,y,r) = y(r)*y(r)*x(r)*x(r)*w(r)
        z_lms[threadId](x, y) += xSquare.transpose() * ySquare;
      } // for itY
    }   // for itX
  }     // for block

  dois = Eigen::MatrixXd::Zero(nX, nY);
#ifdef _OPENMP
  Eigen::setNbThreads(nThreads);
  // sum over all threads
  for (unsigned int t = 0; t < z_lms.size(); ++t) {
    dois += z_lms[t];
  }
#else
  dois = z_lms[0];
#endif
  dois = dois.array().sqrt();
  timeTaken(2, "DOI -- calculations");
}

} /* namespace Serenity */
