/**
 * @file ABCoreHamiltonian.h
 *
 * @date May 15, 2018
 * @author Moritz Bensberg
 *
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOREHAMILTONIAN_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOREHAMILTONIAN_H_

/* Include Serenity Internal Headers */
#include "notification/ObjectSensitiveClass.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/* Forward Declarations */
class Geometry;
class Libint;

/**
 * @class ABCoreHamiltonian ABCoreHamiltonian.h
 * @brief Calculates the outer diagonal block of the core hamiltonian for two basis sets A and B.
 *
 * The matrix entries will have the form\n
 * \f$ h^{AB}_{\nu\mu}=\langle \chi^A_\nu|-\frac{\nabla ^2}{2}-\sum_I \frac{Z_I}{|R_I-r|}| \chi^B_\mu \rangle \f$ ,\n
 * where \f$ \chi^A_\nu \f$ and \f$ \chi^B_\mu \f$ are basis functions of the systems A
 * and B, respectively.
 *
 */
template<Options::SCF_MODES SCFMode>
class ABCoreHamiltonian : public ABPotential<SCFMode>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param geometry The geometry.
   */
  ABCoreHamiltonian(std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB,
                    std::shared_ptr<Geometry> geometry);
  /**
   * @brief Default destructor.
   */
  ~ABCoreHamiltonian() = default;

  /**
   * @brief Getter for the AB fock matrix contribution.
   * @return The AB fock matrix contribution.
   */
  SPMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Deletes the AB fock matrix contribution if it is out of date.
   */
  void notify() override final {
    _abPotential = nullptr;
  };

 private:
  ///@brief The geometry.
  std::shared_ptr<Geometry> _geom;
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief A Libint instance.
  const std::shared_ptr<Libint> _libint;
  ///@brief The ECP contribution.
  std::shared_ptr<ABPotential<SCFMode>> _abEffectiveCorePotential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABCOREHAMILTONIAN_H_ */
