/**
 * @file   MOCalculator.cpp
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#include "data/grid/MOCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "grid/GridController.h"
#include "math/Matrix.h"
#include "misc/HelperFunctions.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

MOCalculator::MOCalculator(std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController)
  : _basisFuncOnGridController(basisFuncOnGridController) {
}

Eigen::MatrixXd MOCalculator::calcMOValuesOnGrid(const Eigen::MatrixXd& coefficients, double mnpTruncationThreshold) {
  takeTime("MO evaluation on grid");
  unsigned int nBlocks = _basisFuncOnGridController->getNBlocks();
  unsigned int nGridPts = _basisFuncOnGridController->getNGridPoints();

  // prepare return value
  const unsigned int nMOs = coefficients.cols();
  Eigen::MatrixXd mosOnGrid = Eigen::MatrixXd::Zero(nGridPts, nMOs);
#ifdef _OPENMP
  Eigen::setNbThreads(1);
#endif
#pragma omp parallel for schedule(dynamic)
  for (unsigned int block = 0; block < nBlocks; block++) {
    const auto& thisBlockData = _basisFuncOnGridController->getBlockOnGridData(block);
    const unsigned int blockFirstIndex = _basisFuncOnGridController->getFirstIndexOfBlock(block);
    const unsigned int blockSize = thisBlockData->functionValues.rows();
    // Basis function value and coefficient prescreening.
    const Eigen::SparseMatrix<double> projection = constructProjectionMatrix(thisBlockData->negligible);
    const Eigen::MatrixXd sigCoefficients_block = projection.transpose() * coefficients;
    const Eigen::MatrixXd sigFunctionValues_block = thisBlockData->functionValues * projection;
    for (unsigned int iMO = 0; iMO < nMOs; ++iMO) {
      const Eigen::VectorXd coefficients_i = sigCoefficients_block.col(iMO);
      const Eigen::SparseMatrix<double> cProjection =
          constructSignificantMnPProjectionMatrix(coefficients_i, mnpTruncationThreshold);
      const Eigen::MatrixXd signCFunctionValues = sigFunctionValues_block * cProjection;
      const Eigen::MatrixXd signCCoefficients = cProjection.transpose() * coefficients_i;
      mosOnGrid.col(iMO).segment(blockFirstIndex, blockSize) = signCFunctionValues * signCCoefficients;
    } // for iMO
  }   // for block
#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif
  timeTaken(3, "MO evaluation on grid");
  return mosOnGrid;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Matrix<double>>
MOCalculator::calcOccMOValuesOnGrid(CoefficientMatrix<SCFMode>& coefficientMatrix,
                                    SpinPolarizedData<SCFMode, unsigned int>& nOccOrbs, double mnpTruncationThreshold) {
  // get some useful values
  unsigned int nGridPts = _basisFuncOnGridController->getNGridPoints();

  // prepare return value
  SpinPolarizedData<SCFMode, Matrix<double>> mosOnGrid;

  for_spin(mosOnGrid, nOccOrbs) {
    mosOnGrid_spin.resize(nGridPts, nOccOrbs_spin);
    mosOnGrid_spin.setZero();
  };

  // call the function for every occupied MO
  for_spin(mosOnGrid, nOccOrbs, coefficientMatrix) {
    const Eigen::MatrixXd coeffs = coefficientMatrix_spin.leftCols(nOccOrbs_spin);
    mosOnGrid_spin = calcMOValuesOnGrid(coeffs, mnpTruncationThreshold);
  };

  return mosOnGrid;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Matrix<double>>
MOCalculator::calcAllMOValuesOnGrid(CoefficientMatrix<SCFMode>& coefficientMatrix, double mnpTruncationThreshold) {
  // call the function above with nBasisFunc as nOccOrbs -> all MOs are returned
  SpinPolarizedData<SCFMode, unsigned int> orbs(_basisFuncOnGridController->getNBasisFunctions());
  SpinPolarizedData<SCFMode, Matrix<double>> mosOnGrid =
      calcOccMOValuesOnGrid<SCFMode>(coefficientMatrix, orbs, mnpTruncationThreshold);

  return mosOnGrid;
}

template SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Matrix<double>>
MOCalculator::calcOccMOValuesOnGrid<Options::SCF_MODES::RESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>& nOccOrbs, double mnpTruncationThreshold);
template SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Matrix<double>>
MOCalculator::calcOccMOValuesOnGrid<Options::SCF_MODES::UNRESTRICTED>(
    CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>& coefficientMatrix,
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>& nOccOrbs, double mnpTruncationThreshold);
template SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Matrix<double>>
MOCalculator::calcAllMOValuesOnGrid<Options::SCF_MODES::RESTRICTED>(CoefficientMatrix<Options::SCF_MODES::RESTRICTED>& coefficientMatrix,
                                                                    double mnpTruncationThreshold);
template SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Matrix<double>>
MOCalculator::calcAllMOValuesOnGrid<Options::SCF_MODES::UNRESTRICTED>(CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>& coefficientMatrix,
                                                                      double mnpTruncationThreshold);

} /* namespace Serenity */
