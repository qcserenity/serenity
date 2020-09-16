/**
 * @file CholeskyDecomposer.h
 *
 * @date Jan 25, 2017
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

#ifndef MATH_LINEARALGEBRA_CHOLESKYDECOMPOSER_H_
#define MATH_LINEARALGEBRA_CHOLESKYDECOMPOSER_H_

/* Include Std and External Headers */
#include <Eigen/Sparse>
#include <functional>
#include <memory>

namespace Serenity {
/**
 * @brief This class holds routines to perform a direct pivoted Cholesky decomposition,
 *        i.e. without ever storing the matrix which is decomposed. This implementation
 *        follows the algorithm of F. Aquilante et. al. To my knowledge, it is the
 *        only direct pivoted Cholesky decomposition that has been published so far.
 *        Most variable names are chosen according to the given reference.
 *
 *        Ref: F. Aquilante, L. Boman, J. Bostroem, H. Koch, R. Lindh, A. S. de Meras
 *             and T. B. Pedersen in Linear-Scaling Techniques in Computational Chemistry
 *             and Physics. Dodrecht: Springer, 2011.
 */
class CholeskyDecomposer {
 public:
  /**
   *
   * @param diagonal
   * @param matrixColumnCalculator
   * @param decompositionThreshold
   * @param screeningDamping
   * @param spanFactor
   * @param nMaxQual
   * @param negativeThreshold
   * @param negativeFailThreshold
   */
  CholeskyDecomposer(Eigen::VectorXd& diagonal,
                     std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<unsigned int>& qualifiedSet)> matrixColumnCalculator,
                     const double decompositionThreshold = 1.0e-8, const double screeningDamping = 1.0,
                     const double spanFactor = 0.01, const unsigned int nMaxQual = 100,
                     const double negativeThreshold = -1.0e-13, const double negativeFailThreshold = -1.0e-8,
                     const double pruningThreshold = 1.0e-10);
  virtual ~CholeskyDecomposer() = default;

  std::unique_ptr<Eigen::SparseMatrix<double>> getCholeskyVectors();

 private:
  // The diagonal elements of the matrix to be decomposed
  Eigen::VectorXd _diagonal;
  // A function calculating qualified columns of the matrix to be decomposed.
  std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<unsigned int>& qualifiedSet)> _matrixColumnCalculator;
  // The decomposition threshold
  const double _decompositionThreshold;
  // The screening damping factor.
  const double _screeningDamping;
  // The span factor
  const double _spanFactor;
  // The maximum number of qualified sets
  const unsigned int _nMaxQual;
  // Diagonal elements which fall below this threshold are set to zero
  const double _negativeThreshold;
  // If diagonal elements fall below this threshold the matrix to be decomposed is
  // probably not positive definite
  const double _negativeFailThreshold;
  // Entries in the Cholesky vectors below this threshold are set to zero (i.e. not stored)
  const double _pruningThreshold;
  // True if decomposition has been performed
  bool _hasBeenCalculated;

  // The Cholesky vectors
  std::unique_ptr<Eigen::SparseMatrix<double>> _choleskyVectors;

  std::vector<unsigned int> calculateQualifiedSet(double dMin);

  void findQmax(double& qMax, unsigned int& qMaxIndex, std::vector<unsigned int>& qualifiedSet);

  void decompose();
};

} /* namespace Serenity */

#endif /* MATH_LINEARALGEBRA_CHOLESKYDECOMPOSER_H_ */
