/**
 * @file DirectCD.h
 *
 * @date Jun 1, 2016
 * @author Michael Boeckers
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

#ifndef SRC_BASICS_MATH_DIRECTCD_H_
#define SRC_BASICS_MATH_DIRECTCD_H_

/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <Eigen/Sparse>

namespace Serenity {
using namespace Eigen;
/**
 * @brief This class holds routines to perform a direct pivoted Cholesky decomposition,
 *        i.e. without ever storing the matrix which is decomposed. This implementation
 *        follows the algorithm of F. Aquilante et. al. To my knowledge, it is the
 *        only direct pivoted Cholesky decomposition that has been published so far.
 *
 *        Ref: F. Aquilante, L. Boman, J. Bostroem, H. Koch, R. Lindh, A. S. de Meras
 *             and T. B. Pedersen in Linear-Scaling Techniques in Computational Chemistry
 *             and Physics. Dodrecht: Springer, 2011.
 */
class DirectCD {
 public:
  /**
   *
   * @param diagonal                The diagonal elements of the matrix
   * @param matrixColumnCalculator  A function which calculates columns of the matrix corresponding
   *                                to a qualified set. The order of the qualified set must not be
   *                                changed.
   * @param decompositionThreshold  The decomposition threshold. The default is 1.0e-8.
   * @param spanFactor              The span factor. This factor modifies the pivoting by allowing
   *                                diagonals to be treated even if the value of the diagonal element
   *                                is not the largest. The default is 0.01.
   */
  DirectCD(VectorXd diagonal, std::function<Matrix<double>(VectorXi& qualifiedSet)> matrixColumnCalculator,
           double decompositionThreshold = 1.0e-8, double spanFactor = 0.01);

  virtual ~DirectCD() = default;
  /**
   * @brief Perform the decomposition
   */
  void decompose();
  /**
   *
   * @return Returns a reference to the Cholesky vectors stored in a sparse matrix
   */
  SparseMatrix<double>& getCholeskyVectors();

 private:
  /**
   * @brief The diagonal elements of the matrix
   */
  VectorXd _diagonal;
  /**
   * @brief A function calculating qualified columns of the matrix to be decomposed.
   */
  std::function<Matrix<double>(VectorXi& qualifiedSet)> _matrixColumnCalculator;
  /**
   * @brief The decomposition threshold
   */
  double _decompositionThreshold = 1.0e-8;
  /**
   * @brief Span factor: The span factor modifies the pivoting by allowing
   *        diagonals to be treated even if the value of the diagonal element
   *        is not the largest.
   */
  double _spanFactor = 0.01;
  /**
   * @brief The Cholesky vectors
   */
  SparseMatrix<double> _choleskyVectors;
  /**
   * @brief The dimension of the matrix to be decomposed
   */
  unsigned int _nDimension;
  /**
   * @brief Pivoting indices
   */
  VectorXi _pivotingIndices;
  /**
   * @brief The number of Cholesky vectors
   */
  unsigned int _nCholesky = 0;
  /**
   * @brief Diagonal elements, which fall below this threshold after each update are set to zero
   */
  double _negativeThreshold = -1.0e-12;
  /**
   * @brief Stop decomposition if diagonal elements become too negative
   */
  double _negativeFailThreshold = -1.0e-8;
  /**
   * @brief               Finds the largest diagonal in qualified set of diagonal elements
   * @param qMax          Largest diagonal element
   * @param qMaxIndex     Index of largest diagonal element
   * @param qualifiedSet  Vector with indices corresponding to qualified diagonals
   */
  void findQmax(double& qMax, unsigned int& qMaxIndex, VectorXi& qualifiedSet);
  /**
   * @brief               Computes a list of qualified diagonals to be decomposed
   * @param dMin          Smallest diagonal element that shall be treated
   * @param qualifiedSet  Qualified set of diagonals
   * @param nQualified    Number of qualified diagonal elements
   */
  void computeQualified(double& dMin, VectorXi& qualifiedSet, unsigned int& nQualified);
};

} /* namespace Serenity */

#endif /* SRC_BASICS_MATH_DIRECTCD_H_ */
