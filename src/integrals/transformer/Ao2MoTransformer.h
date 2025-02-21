/**
 * @file   Ao2MoTransformer.h
 * @author Jan Unsleber
 * @date   23. November 2015
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
#ifndef AO2MOTRANSFORMER_H_
#define AO2MOTRANSFORMER_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "math/Matrix.h"
#include "math/RegularRankFourTensor.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class BasisController;
/**
 * @class Ao2MoTransformer Ao2MoTransformer.h
 * @brief In quantum chemistry, in order to perform calculations, we usually expand the
 *        molecular orbitals (MOs) \f$ \left\{ \phi \right\} \f$ in a finite set of atomic orbitals (AOs)
 *        \f$ \left\{ \chi \right\} \f$
 *        \f[
 *          \phi_i = \sum_\mu c_{i\mu} \chi_\mu \;.
 *        \f]
 *        Most quantities are however needed in the molecular orbital basis and one thus has
 *        to transform from AO to MO basis. While this is easy and computationally cheap for
 *        some quantities, e.g. the one-electron part of the Fock- or Kohn-Sham matrix, it
 *        becomes the time-determining step in post-SCF methods like MP2.
 *        This class provides efficient objects to transform from AO to MO basis .
 */
class Ao2MoTransformer {
 public:
  Ao2MoTransformer(std::shared_ptr<BasisController> basisController);
  virtual ~Ao2MoTransformer() = default;

  /**
   *
   * @param twoElectronIntegrals A rank four tensor holding the two-electron integrals.
   * @param result               A rank four tensor of size nOrbitals which will be filled
   *                             with the transformed integrals. This can be the input
   *                             tensor "twoElectronIntegrals".
   * @param coefficients         The CoefficientMatrix holding the MO coefficients.
   * @param nOrbitals            The number of orbitals to be transformed.
   * @brief                      Transforms all two electron integrals from AO-basis to
   *                             MO-basis.
   *                             The two-electron integral \f$ \left(i j | k l \right) \f$
   *                             can be obtained by
   *                             \f[
   *                            \left(i j | k l \right) = \sum_\mu \sum_\nu \sum_\lambda \sum_\sigma
   *                                                      c_{\mu i} c_{\nu j} c_{\lambda k} c_{\sigma l}
   *                                                      \left(\mu \nu | \lambda \sigma \right) \;.
   *                             \f]
   *                             This expression scales with \f$ N^8 \f$, but can be rewritten as
   *                             \f[
   *                             \left(\mu \nu | \lambda l \right) = \sum_\sigma  c_{\sigma l} \left(\mu \nu | \lambda
   * \sigma \right) \f] \f[ \left(\mu \nu | k l \right) = \sum_\lambda  c_{\lambda k} \left(\mu \nu | \lambda l \right)
   *                             \f]
   *                             \f[
   *                             \left(\mu j | k l \right) = \sum_\nu  c_{\nu j} \left(\mu \nu | k l \right)
   *                             \f]
   *                             \f[
   *                             \left(i j | k l \right) = \sum_\mu  c_{\mu k} \left(\mu j| k l \right) \; ,
   *                             \f]
   *                             where each sum scales with \f$N^5 \f$.
   *                             The transformation is parallelized over the first loop, respectively.
   */
  void transformTwoElectronIntegrals(RegularRankFourTensor<double>& twoElectronIntegrals, RegularRankFourTensor<double>& result,
                                     const Eigen::MatrixXd coefficients, const unsigned int& nOrbitals);

 private:
  unsigned int _nBasisFunc;
};

} // namespace Serenity
#endif /* AO2MOTRANSFORMER_H_ */
