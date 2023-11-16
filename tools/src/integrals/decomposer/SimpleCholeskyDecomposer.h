/**
 * @file   SimpleCholeskyDecomposer.h
 *
 * @date   Mar 2, 2021
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

#ifndef SRC_INTEGRALS_DECOMPOSER_SIMPLECHOLESKYDECOMPOSER_H_
#define SRC_INTEGRALS_DECOMPOSER_SIMPLECHOLESKYDECOMPOSER_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>

namespace Serenity {

/**
 * @class SimpleCholeskyDecomposer SimpleCholeskyDecomposer.h
 *
 *  @brief This class performs a very rudimentary pivoted Cholesky decomposition.
 *  For the algorithm check the manual.
 */
class SimpleCholeskyDecomposer {
 public:
  /**
   * @brief Constructor
   * @param mat The positive semi-definite matrix to decompose.
   * @param thresh The decomposition threshold.
   */
  SimpleCholeskyDecomposer(Eigen::MatrixXd mat, double thresh);

  /**
   *  @brief Getter for the Cholesky vectors
   * @return A Matrix where each column is a Cholesky Vector.
   */
  Eigen::MatrixXd getVectors();

  /**
   *  @brief Getter for the Cholesky basis.
   * @return A vector where the values correspond to the indices of the original matrix
   *                   and the position in the vector corresponds to the Cholesky vector.
   */
  std::vector<unsigned int> getCholeskyBasis();

 private:
  /**
   * @brief Runs the Cholesky decomposition
   * @return True if the decomposition was successful.
   */
  bool decompose();

  double _thresh;
  Eigen::MatrixXd _mat;
  Eigen::MatrixXd _vec;
  Eigen::VectorXd _diag;
  std::vector<unsigned int> _cbas;
};

} /* namespace Serenity */

#endif /* SRC_INTEGRALS_DECOMPOSER_SIMPLECHOLESKYDECOMPOSER_H_ */
