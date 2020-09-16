/**
 * @file   MOCalculator.h
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
#ifndef MOCALCULATOR_H_
#define MOCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/grid/GridData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
template<class T>
class Matrix;
class BasisFunctionOnGridController;
/**
 * @class  Serenity::MOCalculator MOCalculator.h
 * @brief  Calculates the real-space representation of the MOs and the
 *         kinetic energy density on a grid.
 */

class MOCalculator {
 public:
  /**
   * @brief Constructor
   * @param basisFuncOnGridController The BasisFunctionOnGridController providing
   *        the AOs on a grid
   */
  MOCalculator(std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController);
  virtual ~MOCalculator() = default;

  /**
   * @brief Calculates the real-space representation of one MO from the real-space
   *        representation of the AOs.
   * @param coefficients The linear coefficients for the MO to be calculated.
   * @return An Eigen vector containing the grid representations of the MO
   */
  Eigen::MatrixXd calcMOValuesOnGrid(const Eigen::MatrixXd& coefficients, double mnpTruncationThreshold = 1e-5);

  /**
   * @brief Calculates the real-space representation of the occupied MOs using calcMOValuesOnGrid.
   * @param coefficientMatrix The complete CoefficientMatrix containing the
   *        occupied MOs.
   * @param nOccOrbs The number of occupied orbitals.
   * @return A Matrix containing the real-space representation of the occupied MOs
   *         as column vectors.
   */
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, Matrix<double>> calcOccMOValuesOnGrid(CoefficientMatrix<SCFMode>& coefficientMatrix,
                                                                   SpinPolarizedData<SCFMode, unsigned int>& nOccOrbs,
                                                                   double mnpTruncationThreshold = 1e-5);

  /**
   * @brief Calculates the real-space representation of all MOs using calcMOValuesOnGrid.
   * @param coefficientMatrix See above
   * @return A Matrix containing the real-space representation of all MOs
   *         as column vectors.
   */
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, Matrix<double>> calcAllMOValuesOnGrid(CoefficientMatrix<SCFMode>& coefficientMatrix,
                                                                   double mnpTruncationThreshold = 1e-5);

 private:
  std::shared_ptr<BasisFunctionOnGridController> _basisFuncOnGridController;
};

} /* namespace Serenity */

#endif /* MOCALCULATOR_H_ */
