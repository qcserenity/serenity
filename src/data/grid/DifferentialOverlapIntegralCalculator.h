/**
 * @file DifferentialOverlapIntegralCalculator.h
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

#ifndef DATA_GRID_DIFFERENTIALOVERLAPINTEGRALCALCULATOR_H_
#define DATA_GRID_DIFFERENTIALOVERLAPINTEGRALCALCULATOR_H_
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices.
#include <Eigen/SparseCore> //Sparse matrices.
#include <memory>           //smrt_ptr.

namespace Serenity {
/* Forward Declarations */
class BasisFunctionOnGridController;
class ShellPairData;
/**
 * @class DifferentialOverlapIntegralCalculator DifferentialOverlapIntegralCalculator.h
 * @brief Calculates the differential overlap between two sets of functions \f$ \{X\},\{Y\} \f$ expressed
 *        over basis functions.\n
 *        \f$ DOI(X,Y) = \sqrt(\int X(r)^2 Y(r)^2 \mathrm{d}r \f$,\n
 *        where \f$ X \f$ and \f$ Y \f$ are given by\n
 *        \f$ X(r) = \sum_{\nu} c_{x\nu} \chi_\nu \f$\n
 *        and \f$ \chi_\nu \f$ is a basis function with its respective linear coefficients \f$ c_{x\nu} \f$.\n\n
 *
 *        This implementation is linear scaling.
 */
class DifferentialOverlapIntegralCalculator {
 private:
  /**
   * Purely static never instantiated!
   */
  DifferentialOverlapIntegralCalculator() = default;
  virtual ~DifferentialOverlapIntegralCalculator() = default;

 public:
  /**
   * @brief Calculates the differential overlap integrals (DOIs) between function expressed
   *        as linear combinations of basis functions.
   * @param c_x Linear coefficients of the first function.
   * @param c_y Linear coefficients of the second function.
   * @param basisFunctionToXMap Prescreening map for X.
   * @param basisFunctionToYMap Prescreening map for Y.
   * @param basisFuncOnGridController The basis function on grid controller.
   * @param dois Results are written into this object. DOIs are sorted as X-->row and Y-->column.
   */
  static void calculateDOI(const Eigen::MatrixXd& c_x, const Eigen::MatrixXd& c_y,
                           const Eigen::SparseMatrix<int>& basisFunctionToXMap,
                           const Eigen::SparseMatrix<int>& basisFunctionToYMap,
                           std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController, Eigen::MatrixXd& dois);

  /**
   * @brief Calculates the differential overlap integrals (DOIs) between basis functions on a grid.
   * @param basisFuncOnGridController The basis function on grid controller.
   * @param dois Results are written into this object. DOIs are sorted as X-->row and Y-->column.
   */
  static void calculateDOI(std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController, Eigen::MatrixXd& dois);

  /**
   * @brief Generate the shell pair data object for integral prescreening based on DOI values.
   * @param basisFuncOnGridController The basis function on grid controller.
   * @param cutOff The integral cut off. By default, 1e-12.
   * @return The shell pair data objects.
   */
  static std::vector<ShellPairData>
  calculateDOIShellPairData(std::shared_ptr<BasisFunctionOnGridController> basisFuncOnGridController, double cutOff = 1e-12);
};

} /* namespace Serenity */

#endif /* DATA_GRID_DIFFERENTIALOVERLAPINTEGRALCALCULATOR_H_ */
