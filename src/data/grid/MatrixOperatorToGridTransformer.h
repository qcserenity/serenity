/**
 * @file MatrixOperatorToGridTransformer.h
 *
 * @date May 15, 2015
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
#ifndef MATRIXOPERATORTOGRIDTRANSFORMER_H
#define MATRIXOPERATORTOGRIDTRANSFORMER_H
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGrid.h"
#include "data/grid/GridData.h"
#include "data/matrices/MatrixInBasis.h"
#include "math/Derivatives.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <functional>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class BasisFunctionOnGridController;

/**
 * @class MatrixOperatorToGridTransformer MatrixOperatorToGridTransformer.h
 * @brief Holds methods to discretize an operator in matrix representation on grid points
 *
 * Caution!! Using this class only makes sense for certain matrices, because the contributions
 * of the different entries in the matrix are effectively added up. This is correct to e.g.
 * transform the DensityMatrix to an integration grid; but it is wrong for the transformation
 * of a potential or the FockMatrix! Check carefully, whether the formula given in the method
 * description is actually what you want.
 *
 * Purely static class.
 *
 * For tests: The DensityOnGridCalculator does not much more than using the routines
 * provided here, so take a look at the tests of that class.
 */
class MatrixOperatorToGridTransformer {
 public:
  MatrixOperatorToGridTransformer() = delete;
  /**
   * @brief Calculates \f$ {\rm result}(r) = \sum_{i,j} \chi_i(r) \cdot {\rm matrix}_{i,j} \chi_j(r) \f$ for all grid
   * points.
   *
   * @param matrix representation of the operator;
   *    not every operator is suited to be used with this function!
   * @param result the discretized values on the grid (see the formula above)
   * @param resultGradient if specified it holds the gradient of result on exit:\n
   *        \f$ {\rm resultGradient}(r) =
   *               \sum_{i,j} \chi_i(r) \cdot {\rm matrix}_{i,j} \nabla\chi_j(r) +
   *               \sum_{i,j} \nabla\chi_i(r) \cdot {\rm matrix}_{i,j} \chi_j(r) \f$ for all grid points
   * @param resultHessian  if specified it holds the Hessian of result on exit:\n
   *        \f$ {\rm resultHessian}_{a,b}(r) =
   *               \sum_{i,j} \chi_i(r) \cdot {\rm matrix}_{i,j}
   *                          \frac{\partial}{\partial a}\frac{\partial}{\partial b}\chi_j(r) +
   *               \sum_{i,j} \frac{\partial}{\partial a}\frac{\partial}{\partial b}\chi_i(r) \cdot
   *                           {\rm matrix}_{i,j} \chi_j(r) +
   *               \sum_{i,j} \frac{\partial}{\partial a}\chi_i(r) \cdot {\rm matrix}_{i,j}
   *                          \frac{\partial}{\partial b}\chi_j(r) \f$ for all grid points +
   *               \sum_{i,j} \frac{\partial}{\partial b}\chi_i(r) \cdot {\rm matrix}_{i,j}
   *                          \frac{\partial}{\partial a}\chi_j(r) \f$ for all grid points \f$
   * @param basisFunctionOnGridController working on the basis in which matrix is defined and giving
   *    data on the grid for which result shall be calculated
   * @param blockAverageThreshold if the average contribution of a block of grid points is below
   *                              this threshold it is skipped for efficiency
   * Caution: The derivative operator is only added to the basis functions. No fancy magic
   * happens to the data in matrix.\n
   *
   * Several overloads for this function exist to make the calculation of the Hessian (or of gradient
   * and Hessian) optional, and to treat several matrices in one run (for computational efficiency).
   *
   */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<RESTRICTED>& matrix, Eigen::VectorXd& result,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrices = {matrix};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> results = {result};
    return transform(matrices, basisFunctionOnGridController, results);
  }
  /** See overload above, single matrix, values and gradient */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<RESTRICTED>& matrix, DensityOnGrid<RESTRICTED>& result,
                                            Gradient<DensityOnGrid<RESTRICTED>>& resultGradient,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrices = {matrix};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> results = {result};
    std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>> resultGradients = {
        {resultGradient.x, resultGradient.y, resultGradient.z}};
    return transform(matrices, basisFunctionOnGridController, results, &resultGradients);
  }
  /** See overload above, single matrix, values, gradient and Hessian */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<RESTRICTED>& matrix, DensityOnGrid<RESTRICTED>& result,
                                            Gradient<DensityOnGrid<RESTRICTED>>& resultGradient,
                                            Hessian<DensityOnGrid<RESTRICTED>>& resultHessian,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrices = {matrix};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> results = {result};
    std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>> resultGradients = {
        {resultGradient.x, resultGradient.y, resultGradient.z}};
    std::vector<Hessian<std::reference_wrapper<Eigen::VectorXd>>> resultHessians = {
        {resultHessian.xx, resultHessian.xy, resultHessian.xz, resultHessian.yy, resultHessian.yz, resultHessian.zz}};
    return transform(matrices, basisFunctionOnGridController, results, &resultGradients, &resultHessians);
  }
  /**
   * @brief see the other method. For the same basis and the same grid two matrices are transformed.
   *
   * This saves some calculation time in contrast to calling the other method twice. Used e.g. for
   * spin-polarized data.
   *
   * @param matrices
   * @param results
   * @param basisFunctionOnGridController
   * @param derivativeOrder
   *
   */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<UNRESTRICTED>& matrices, DensityOnGrid<UNRESTRICTED>& results,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrixVec = {matrices.alpha, matrices.beta};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> resultVec = {results.alpha, results.beta};
    return transform(matrixVec, basisFunctionOnGridController, resultVec);
  }
  /** See overload above, two matrices, values and gradients */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<UNRESTRICTED>& matrices, DensityOnGrid<UNRESTRICTED>& results,
                                            Gradient<DensityOnGrid<UNRESTRICTED>>& resultGradients,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrixVec = {matrices.alpha, matrices.beta};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> resultVec = {results.alpha, results.beta};
    std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>> resultGradientVec = {
        {resultGradients.x.alpha, resultGradients.y.alpha, resultGradients.z.alpha},
        {resultGradients.x.beta, resultGradients.y.beta, resultGradients.z.beta}};
    return transform(matrixVec, basisFunctionOnGridController, resultVec, &resultGradientVec);
  }
  /** See overload above, two matrices, values, gradients and Hessians */
  static Eigen::SparseVector<int> transform(const MatrixInBasis<UNRESTRICTED>& matrices, DensityOnGrid<UNRESTRICTED>& results,
                                            Gradient<DensityOnGrid<UNRESTRICTED>>& resultGradients,
                                            Hessian<DensityOnGrid<UNRESTRICTED>>& resultHessians,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> matrixVec = {matrices.alpha, matrices.beta};
    std::vector<std::reference_wrapper<Eigen::VectorXd>> resultVec = {results.alpha, results.beta};
    std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>> resultGradientVec = {
        {resultGradients.x.alpha, resultGradients.y.alpha, resultGradients.z.alpha},
        {resultGradients.x.beta, resultGradients.y.beta, resultGradients.z.beta}};
    std::vector<Hessian<std::reference_wrapper<Eigen::VectorXd>>> resultHessianVec = {
        {resultHessians.xx.alpha, resultHessians.xy.alpha, resultHessians.xz.alpha, resultHessians.yy.alpha,
         resultHessians.yz.alpha, resultHessians.zz.alpha},
        {resultHessians.xx.beta, resultHessians.xy.beta, resultHessians.xz.beta, resultHessians.yy.beta,
         resultHessians.yz.beta, resultHessians.zz.beta}};
    return transform(matrixVec, basisFunctionOnGridController, resultVec, &resultGradientVec, &resultHessianVec);
  }
  /**
   * @brief see the other methods. For the same basis and the same grid many matrices are transformed.
   *
   * This saves some calculation time in contrast to calling the other method several times.
   *
   * @param matrices
   * @param results
   * @param basisFunctionOnGridController
   * @param derivativeOrder
   */
  static Eigen::SparseVector<int> transform(const std::vector<std::reference_wrapper<const MatrixInBasis<RESTRICTED>>>& matrices,
                                            std::vector<std::reference_wrapper<Eigen::VectorXd>>& results,
                                            BasisFunctionOnGridController& basisFunctionOnGridController) {
    std::vector<std::reference_wrapper<const Eigen::MatrixXd>> symMatrices;
    for (const auto& matrix : matrices)
      symMatrices.push_back(matrix.get());
    std::vector<std::reference_wrapper<Eigen::VectorXd>> vectors;
    for (const auto& result : results)
      vectors.push_back(result.get());
    return transform(symMatrices, basisFunctionOnGridController, vectors);
  }

 private:
  // Generalized private method. This actually performs the computations and is called by the others.
  static Eigen::SparseVector<int>
  transform(const std::vector<std::reference_wrapper<const Eigen::MatrixXd>>& matrices,
            BasisFunctionOnGridController& basisFunctionOnGridController,
            std::vector<std::reference_wrapper<Eigen::VectorXd>>& results,
            std::vector<Gradient<std::reference_wrapper<Eigen::VectorXd>>>* resultGradientsPtr = nullptr,
            std::vector<Hessian<std::reference_wrapper<Eigen::VectorXd>>>* resultHessiansPtr = nullptr);
};

} /* namespace Serenity */
#endif /* MATRIXOPERATORTOGRIDTRANSFORMER_H */
