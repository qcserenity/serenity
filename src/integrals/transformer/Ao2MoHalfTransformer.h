/**
 * @file   Ao2MoHalfTransformer.h
 * @author Moritz Bensberg
 * @date   17. December 2018
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
#ifndef AO2MOHALFTRANSFORMER_H_
#define AO2MOHALFTRANSFORMER_H_

/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices.
#include <memory>      //smrt_ptr.

namespace Serenity {
/* Forward declarations */
class BasisController;
/**
  @class Ao2MoHalfTransformer Ao2MoHalfTransformer..h
  @brief Calculates the integrals (ri|sj) for a given i and j. The indices i and j
         denote occupied orbitals and the indices r and s denote AO basis functions.
         \n\n

         \f$ (ri|sj) = \sum_{\nu\mu} c_{i\mu}c_{j\nu}(r\mu|s\nu) \f$ \n
*/
class Ao2MoHalfTransformer {
 public:
  /**
  @brief Constructor.
  @param basisControllerA The basis over which i and j are expanded.
  @param basisControllerB The basis over which r and s is expanded.
  */
  Ao2MoHalfTransformer(std::shared_ptr<BasisController> basisControllerA, std::shared_ptr<BasisController> basisControllerB)
    : _basisControllerA(basisControllerA), _basisControllerB(basisControllerB) {
  }
  virtual ~Ao2MoHalfTransformer() = default;
  /**
   * @brief Calculate Fock like exchange matrix.
   * @param result The integrals.
   * @param pairDensityMatrix The pair-density matrix.
   */
  void transformTwoElectronIntegrals(Eigen::MatrixXd& result, const Eigen::MatrixXd& pairDensityMatrix);

 private:
  // The basis controller A.
  std::shared_ptr<BasisController> _basisControllerA;
  // The basis controller B.
  std::shared_ptr<BasisController> _basisControllerB;
};

} // namespace Serenity
#endif /* AO2MOHALFTRANSFORMER_H_ */
