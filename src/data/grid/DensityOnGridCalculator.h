/**
 * @file   DensityOnGridCalculator.h
 *
 * @date   Mar 15, 2014
 * @author Dennis Barton, severely modified by Thomas Dresselhaus
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
#ifndef DENSITYONGRIDCALCULATOR_H_
#define DENSITYONGRIDCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGrid.h"
#include "data/matrices/DensityMatrix.h"
#include "math/Derivatives.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <iostream>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class BasisFunctionOnGridController;
class GridController;
/**
 * @class DensityOnGridCalculator DensityOnGridCalculator.h
 *
 * @brief Build up electron density on given grid
 *
 * This class provides the density on a precalculated grid
 * depending on the basis and grid of the system.
 */
template<Options::SCF_MODES SCFMode>
class DensityOnGridCalculator {
 public:
  /**
   * @param basisFunctionOnGridController Controller for the basis functions on the grid.
   * @param blockAverageThreshold If the average contribution of a block of grid points is below
   *                              this threshold it is skipped for efficiency.
   */
  DensityOnGridCalculator(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                          double blockAverageThreshold);

  virtual ~DensityOnGridCalculator() = default;
  /**
   * @brief Uses the class MatrixOperatorToGridTransformer to calculate
   * \f$ {\rm result}(r) = \sum_{i,j} \chi_i(r) \cdot {\rm matrix}_{i,j} \chi_j(r) \f$ for all grid
   * points.
   *
   * @param[in] densityMatrix the electron density of the system in matrix form
   * @returns the density on the grid as vector
   */
  DensityOnGrid<SCFMode> calcDensityOnGrid(const DensityMatrix<SCFMode>& densityMatrix);
  /**
   * Call this to avoid a copy of the density on the grid
   * @brief Uses the class MatrixOperatorToGridTransformer to calculate
   * \f$ {\rm result}(r) = \sum_{i,j} \chi_i(r) \cdot {\rm matrix}_{i,j} \chi_j(r) \f$ for all grid
   * points.
   *
   * @param[in]  densityMatrix the electron density of the system in matrix form
   * @param[out] densityOnGrid the density on the grid as vector, will be resized/overwritten
   */
  void calcDensityOnGrid(const DensityMatrix<SCFMode>& densityMatrix, DensityOnGrid<SCFMode>& densityOnGrid);
  /**
   * @brief Uses the class MatrixOperatorToGridTransformer to calculate
   * the density and its gradient for all grid points.
   *
   * @param[in] densityMatrix the electron density of the system in matrix form
   * @param[out] densityGradientOnGrid gradient additional to the density also the density gradient is calculated when
   *                      calling this method.
   * @returns the density on the grid as vector
   */
  DensityOnGrid<SCFMode> calcDensityAndGradientOnGrid(const DensityMatrix<SCFMode>& densityMatrix,
                                                      Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid);
  /**
   * Call this to avoid a copy of the density on the grid
   * @brief Uses the class MatrixOperatorToGridTransformer to calculate
   * the density and its gradient for all grid points.
   *
   * @param[in] densityMatrix the electron density of the system in matrix form
   * @param[out] densityOnGrid the density on the grid as vector, will be resized/overwritten
   * @param[out] densityGradientOnGrid gradient additional to the density also the density gradient is calculated when
   *                      calling this method.
   */
  void calcDensityAndGradientOnGrid(const DensityMatrix<SCFMode>& densityMatrix, DensityOnGrid<SCFMode>& densityOnGrid,
                                    Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid);
  /**
   * @brief Calculates the density and its first and second derivatives on the grid.
   *
   * @param[in] densityMatrix the electron density of the system in matrix form
   * @param[out] densityOnGrid the density on the grid as vector, will be resized/overwritten
   * @param[out] densityGradientOnGrid gradient additional to the density also the density gradient is calculated when
   *                      calling this method.
   * @param[out] densityHessianOnGrid Hessian additional to the density also the density Hessian is calculated when
   *                      calling this method.
   */
  void calcDensityAndDerivativesOnGrid(const DensityMatrix<SCFMode>& densityMatrix, DensityOnGrid<SCFMode>& densityOnGrid,
                                       Gradient<DensityOnGrid<SCFMode>>& densityGradientOnGrid,
                                       Hessian<DensityOnGrid<SCFMode>>& densityHessianOnGrid);
  /**
   * @returns the used BasisFunctionOnGridController.
   */
  inline std::shared_ptr<BasisFunctionOnGridController> getBasisFunctionOnGridController() const {
    return _basisFunctionOnGridController;
  }

  /**
   * @returns the Grid(controller) the object is working on
   */
  std::shared_ptr<GridController> getGridController() const;
  /**
   * @brief Getter for the grid blocks with values different from zero.
   */
  Eigen::SparseVector<int> getNonNegligibleBlocks() {
    return _nonNegligible;
  }

 private:
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridController;
  const double _blockAverageThreshold;
  Eigen::SparseVector<int> _nonNegligible;
};

} /* namespace Serenity */

#endif /* DENSITYONGRIDCALCULATOR_H_ */
