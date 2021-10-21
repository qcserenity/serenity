/**
 * @file ABEffectiveCorePotential.h
 *
 * @date Dec 4, 2018
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABEFFECTIVECOREPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABEFFECTIVECOREPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/* Forward Declarations */
class Atom;

/**
 * @class ABEffectiveCorePotential ABEffectiveCorePotential.h
 * @brief An ABPotential that manages the ECP contribution to the off diagonal block.
 */
template<Options::SCF_MODES SCFMode>
class ABEffectiveCorePotential : public ABPotential<SCFMode>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param atoms The atoms.
   */
  ABEffectiveCorePotential(std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB,
                           std::vector<std::shared_ptr<Atom>> atoms);
  /**
   * @brief Default destructor.
   */
  virtual ~ABEffectiveCorePotential() = default;

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
  ///@brief The atoms.
  std::vector<std::shared_ptr<Atom>> _atoms;
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief
  bool _notZero = false;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABEFFECTIVECOREPOTENTIAL_H_ */
