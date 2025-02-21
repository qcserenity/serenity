/**
 * @file BasisFunctionMapper.h
 *
 * @date 10 Aug 2019
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

#ifndef BASIS_BASISFUNCTIONMAPPER_H_
#define BASIS_BASISFUNCTIONMAPPER_H_
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <memory>

namespace Serenity {
/* Forward declarations */
class BasisController;
class Shell;

/**
 * @class BasisFunctionMapper BasisFunctionMapper.h
 * @brief A class that maps between basis functions of different basis sets.
 */
class BasisFunctionMapper {
 public:
  /**
   * @brief Constructor.
   * @param basisControllerA The basis set to which all further sets are mapped to.
   */
  BasisFunctionMapper(std::shared_ptr<BasisController> basisControllerA);
  /**
   * @brief Default destructor.
   */
  ~BasisFunctionMapper() = default;
  /**
   * @brief Getter for a sparse projection/sorting matrix, that allows mapping between basis functions
   *        of different basis sets. The matrix will look like this:\n\n
   *
   * Example for basis B: B B A1 A3 B A4 shells\n
   * Example for basis A: A1 ... A5\n
   *                      A1 A2 A3 A4 A5\n
   *                  B   0  0  0  0  0\n
   *                  B   0  0  0  0  0\n
   *                  A1  1  0  0  0  0\n
   *                  A3  0  0  1  0  0\n
   *                  B   0  0  0  0  0\n
   *                  A4  0  0  0  1  0\n\n
   * This matrix has the property that multiplying it to a matrix of type basisA x something
   * extracts the rows mapped to the basis B. For an example application @see HuzinagaProjectionPotential.cpp
   * @param basisControllerB The basis controller to which A is compared to.
   * @return
   */
  std::shared_ptr<Eigen::SparseMatrix<double>> getSparseProjection(std::shared_ptr<BasisController> basisControllerB);
  /**
   * @brief Getter for the basis set spanned by the shells {B}\\{A}.
   * @param basisControllerB The basis controller to which A is compared to.
   * @return The differential basis of A and B or an nullptr if both basis sets are identical.
   */
  std::shared_ptr<BasisController> getDifferentialBasis(std::shared_ptr<BasisController> basisControllerB);
  /**
   * @brief Getter for the basis set spanned by the shells {A} U {B}
   * @param basisControllerB The basis controller with which A is combined.
   * @return The combined basis of A and B. For identical basis sets this is A.
   */
  std::shared_ptr<BasisController> getCombinedBasis(std::shared_ptr<BasisController> basisControllerB);

 private:
  // Searches for other in the basis A and returns its index. The index will be larger than
  // the total number of shells if other is not part of A.
  unsigned int getShellIndex(const Shell& other);
  // The basis controller A.
  std::shared_ptr<BasisController> _basisControllerA;
};

} /* namespace Serenity */

#endif /* BASIS_BASISFUNCTIONMAPPER_H_ */
