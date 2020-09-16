/**
 * @file   ScalarOperatorToMatrixAdder.h
 *
 * @date   Mar 22, 2014
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
#ifndef SCALAROPERATORTOMATRIXADDER_H_
#define SCALAROPERATORTOMATRIXADDER_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/GridPotential.h"
#include "data/matrices/FockMatrix.h"
#include "math/Derivatives.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <omp.h>
#include <iostream>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class BasisFunctionOnGridController;
/**
 * @class  ScalarOperatorToMatrixAdder ScalarOperatorToMatrixAdder.h
 * @brief  Transforms some data on a grid into a matrix representation for a given basis (grid integration)
 *
 * The numerical integrals solved here are of the form:\n
 * \f$ \left<\chi_{\mu}(r) | op(r) | \chi_{\nu}(r)\right> \approx
 *     \sum \chi_{\mu}(r) op(r) weight(r) \chi_{\nu}(r) \f$, where weight(r) is the integration
 * weight of the corresponding grid point.
 */
template<Options::SCF_MODES SCFMode>
class ScalarOperatorToMatrixAdder {
 public:
  /**
   * @param basisFunctionOnGridController
   * @param blockAveThreshold if the average contribution of a block of grid points is below
   *                          this threshold it is skipped for efficiency
   */
  ScalarOperatorToMatrixAdder(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridController,
                              double blockAveThreshold);

  ScalarOperatorToMatrixAdder(std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerA,
                              std::shared_ptr<BasisFunctionOnGridController> basisFunctionOnGridControllerB,
                              double blockAveThreshold);

  virtual ~ScalarOperatorToMatrixAdder() {
#ifdef _OPENMP
    omp_destroy_lock(&_lock);
#endif
  }
  /**
   * @brief Adds contributions from a scalar operator which is represented on a grid to a matrix.
   *        Generalized version for SPMatrix.
   *
   * @param matrix holds data expressed in combinations of basis functions. E.g. the Fock matrix.
   * @param scalarOperator expressed on a grid. This vector op will contribute with
   *                       \f$ F_{\mu, \nu} += \left<\chi_{\mu}(r) | op(r) | \chi_{\nu}(r)\right> \f$
   */
  void addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix, const GridPotential<SCFMode>& scalarOperator);

  /**
   * @brief Adds contributions from a scalar operator which is represented on a grid to a matrix.
   *        Generalized version for SPMatrix.
   *
   * @param matrix holds data expressed in combinations of basis functions. E.g. the Fock matrix.
   * @param scalarOperator expressed on a grid. This vector op will contribute with
   *                       \f$ F_{\mu, \nu} += \left<\chi_{\mu}(r) | op(r) | \chi_{\nu}(r)\right> \f$
   * @param gradientOperator This vector g will contribute with
   *                         \f$ F_{\mu, \nu} += \left<\chi_{\mu} | g  | \nabla \chi_{\nu}\right> + \left<\nabla
   * \chi_{\mu} | g | \chi_{\nu}\right> \f$, optional.
   */
  void addScalarOperatorToMatrix(SPMatrix<SCFMode>& matrix, const GridPotential<SCFMode>& scalarOperator,
                                 const Gradient<GridPotential<SCFMode>>& gradientOperator);
  /**
   * @returns the used BasisFunctionOnGridController; Determines Grid and Basis.
   */
  inline std::shared_ptr<BasisFunctionOnGridController> getBasisFunctionOnGridController() const {
    return _basisFunctionOnGridControllerA;
  }

 private:
  /// @brief Basis function value controller for basis A.
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridControllerA;
  /// @brief Basis function value controller for basis B.
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridControllerB;
  /// @brief Block average threshold for preescreening of average basis function values in a block.
  const double _blockAveThreshold;

  /**
   * @brief Integrates over one block.
   * @param iBlock The index of the block.
   * @param blockDataA The basis function values of block A.
   * @param blockDataB The basis function values of block B.
   * @param m_AB The matrix where the result is added to.
   * @param scalarPart The scalar potential for this block.
   */
  void addBlock(unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
                std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB,
                SPMatrix<SCFMode>& m_AB, const GridPotential<SCFMode>& scalarPart);

  /**
   * @brief Integrates over one block.
   * @param iBlock The index of the block.
   * @param blockDataA The basis function values of block A.
   * @param blockDataB The basis function values of block B.
   * @param m_AB The matrix where the result is added to.
   * @param scalarPart The scalar potential for this block.
   * @param gradientPart The gradient part of the potential for this block.
   */
  void addBlock(unsigned int iBlock, std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataA,
                std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockDataB,
                SPMatrix<SCFMode>& m_AB, const GridPotential<SCFMode>& scalarPart,
                const Gradient<GridPotential<SCFMode>>& gradientPart);

#ifdef _OPENMP
  omp_lock_t _lock;
#endif
};

} /* namespace Serenity */

#endif /* SCALAROPERATORTOMATRIXADDER_H_ */
