/**
 * @file   CholeskyDecomposer.h
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef CHOLESKYDECOMPOSER_H_
#define CHOLESKYDECOMPOSER_H_

/* Include Serenity Internal Headers */
#include "integrals/CDIntegralController.h"
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <Eigen/Eigen>

namespace Serenity {

/**
 * @class CholeskyDecomposer CholeskyDecomposer.h
 *
 * @brief The class that performs the pivoted, incomplete Cholesky decomposition of a
 * given matrix. The matrix is provided via the lambda function
 * matrixColumnCalculator, which returns elements of the matrix that correspond
 * to the qualifiedSet(rows) and reducedSet(columns).
 *
 * The pivoting algorithm orders the columns to decompose according to their
 * significance(magnitude of diagonal element) and assures numerical stability.
 *
 * The incomplete decomposition guarantees a fast decomposition combined with a
 * guarantee that the original matrix can be reconstructed with an error that is
 * smaller than the _decompositionThreshold.
 *
 * This corresponds to a modification of the algorithm described in
 * "Francesco Aquilante, T. P., Linear-scaling Techniques in Computational Chem-
 * istry and Physics; Springer Netherlands: 2011; Vol. 13, pp 301–343"
 *
 * For the detailed algorithm check the manual.
 *
 */
class CholeskyDecomposer {
 public:
  /**
   * @brief Constructor.
   *
   * @param matrixName:             A label to identify the matrix.
   * @param cdIntegralController:   The pointer to the cdIntegralController of the used system.
   * @param diagonal:               The diagonal of the matrix written as a vector.
   * @param matrixColumnCalculator: Lambda function that calculates matrix elements. The elements are defined
   *                                                    by the vectors passed to this function (The length of the
   * vectors is the number of columns/rows in the full matrix. If an element is not present in a set the corresponding
   * value in the vector is -1. If the element is present its value corresponds to its position in the set. qualIndices:
   * corresponds to the qualified set or row indices redIndices: corresponds to the reduced set or column indices
   *
   *     // clang-format off
   *  Full Matrix: 1 2 3       redIndices:  1 -1 0        Calculated Block: 3 1
   *               4 5 6       qualIndices: 0 1 2                           6 4
   *               7 8 9                                                    9 7
   *     // clang-format on
   *
   * @param decompositionThreshold: The threshold of the decomposition that determines its accuracy.
   * @param screeningDamping:       Damping factor of the screening in the decomposition.
   * @param spanFactor:             Range of diagonal elements treated in one iteration.
   * @param nMaxQual:               Maximum number of columns that are treated at once.
   * @param negativeThreshold:      Threshold to damp decomposition for semi-definite matrices.
   * @param negativeFailThreshold:  Threshold to fail decomposition if the matrix is not positive (semi-)definite.
   */
  CholeskyDecomposer(std::string matrixName, std::shared_ptr<CDIntegralController> cdIntegralController,
                     Eigen::VectorXd& diagonal,
                     std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<int>& qualIndices, std::vector<int>& redIndices)> matrixColumnCalculator,
                     const double decompositionThreshold = 1.0e-8, const double screeningDamping = 1.0,
                     const double spanFactor = 0.01, const unsigned int nMaxQual = 100,
                     const double negativeThreshold = -1.0e-13, const double negativeFailThreshold = -1.0e-7);

  /// @brief Default destructor.
  virtual ~CholeskyDecomposer() = default;

  /**
   * @brief Performs the Cholesky decomposition
   */
  void run();

  /**
   * @brief getter for the Cholesky Basis. Performs the decomposition if necessary.
   * @return Returns a vector containing a value for each Cholesky vector. The value is the
   *                   the corresponding column index of the original matrix.
   */
  std::vector<unsigned int> getCholeskyBasis();

 private:
  /**
   * @brief Finds the maximum diagonal element in the qualified set
   *
   * @param qMax:         The maximum element in the qualified set.
   * @param qMaxIndex:    The index of the maximum element in the qualified set.
   * @param qualifiedSet: The set of qualified elements.
   */
  void findQmax(double& qMax, unsigned int& qMaxIndex, std::vector<unsigned int>& qualifiedSet);

  /**
   * @brief Performs the Cholesky decomposition of any positive semi-definite
            hermitian matrix
   *
   * @param label: The label to identify the corresponding matrix
   */
  void decompose(std::string label);

  // Label to identify the matrix being decomposed
  std::string _matrixName;
  // decomposition status
  bool _decomposed;
  // The Cholesky integral Controller for the system the matrix belongs to.
  std::shared_ptr<CDIntegralController> _cdIntegralController;
  // The diagonal of the matrix to be decomposed
  Eigen::VectorXd _diagonal;
  // A function calculating qualified columns of the matrix to be decomposed.
  std::function<std::unique_ptr<Eigen::MatrixXd>(std::vector<int>& qualIndices, std::vector<int>& redIndices)> _matrixColumnCalculator;
  // The threshold for the Cholesky decomposition
  const double _decompositionThreshold;
  // The screening damping parameter
  const double _screeningDamping;
  // The factor giving the range of diagonal elements treated at once
  // dMin = _spanFactor * dMax
  const double _spanFactor;
  // The maximum size of the qualified set
  unsigned int _nMaxQual;
  // Diagonal elements which fall below this threshold are set to zero
  const double _negativeThreshold;
  // If diagonal elements fall below this threshold, the matrix to be decomposed is
  // probably not positive definite
  const double _negativeFailThreshold;
  // Current number of calculated Cholesky vectors
  unsigned int _maxJ;
  // The Cholesky Basis used in the Decomposition
  std::vector<unsigned int> _cBasis;
};

} /* namespace Serenity */

#endif /* CHOLESKYDECOMPOSER_H_ */
