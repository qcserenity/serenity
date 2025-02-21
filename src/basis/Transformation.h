/**
 * @file   Transformation.h
 *
 * @date   08-28-2013, 07-08-2016
 * @author Thomas Dresselhaus, Michael Boeckers
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
#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class BasisController;
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;
class OneElectronIntegralCalculator;
template<Options::SCF_MODES SCFMode>
class OrbitalController;

/**
 * @class Transformation Transformation.h
 * @brief Transformation of molecular orbitals defined in a basis A to molecular orbitals
 *        defined in a basis B. This routine is e.g. used for the extended Hueckel guess
 *        to transform the Hueckel orbitals, which are defined in a minimal basis set, into
 *        the basis which is actually used in the calculation.
 *        The transformation can be done by defining a projection operator \f$ P_{BA} \f$
 *        which projects the molecular orbitals in basis A onto basis B. The projection
 *        operator is given by
 *        \f[
 *        \mathbf{P}_{BA} = \mathbf{S}_B^{-1} \mathbf{S}_{AB} \; ,
 *        \f]
 *        where \f$ \mathbf{S}_B^{-1}\f$ is the inverse of the atomic orbital overlap integrals
 *        of basis B and \f$\mathbf{S}_{AB}\f$ is the mixed overlap between atomic orbitals
 *        of basis A and B.
 *        The coefficients of the molecular orbitals in basis B can then be obtained by
 *        \f[
 *        \mathbf{C}_B = \mathbf{P}_{BA} \mathbf{C}_A \; .
 *        \f]
 *
 *        Ref: For the form of the projection operator see
 *             G. Knizia, J. Chem. Theory. Comput. 9, 4834 (2013).
 *             It's derivation can be found in
 *             J. W. Boughton, P. Pulay, J. Comput. Chem. 14, 736 (1993)
 */
class Transformation {
 private:
  /**
   * @brief Default constructor. Never instantiated. Purely static.
   */
  Transformation() = default;

 public:
  /**
   * @brief Default destructor.
   */
  virtual ~Transformation() = default;
  /**
   * @brief Transform molecular orbitals defined in a basis A into molecular orbitals defined in a basis B
   * @param orbitalsA         Molecular orbitals defined in basis A
   * @param basisControllerB  The basis controller of basis B
   * @param overlapB          The overlap integrals of basis B, i.e. \f$\mathbf{S}_B\f$
   * @return                  Returns new molecular orbitals defined in basis B
   */
  template<Options::SCF_MODES SCFMode>
  static std::unique_ptr<OrbitalController<SCFMode>> transformMOs(OrbitalController<SCFMode>& orbitalsA,
                                                                  std::shared_ptr<BasisController> basisControllerB,
                                                                  const MatrixInBasis<RESTRICTED>& overlapB);
};

} /* namespace Serenity */
#endif /* TRANSFORMATION_H_ */
