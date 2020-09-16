/**
 * @file ABZeroPotential.h
 *
 * @date Jun 25, 2018
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABZEROPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABZEROPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/**
 * @class ABZeroPotential ABZeroPotential.h
 * @brief Zero/Dummy ABPotential.
 */
template<Options::SCF_MODES SCFMode>
class ABZeroPotential : public ABPotential<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   */
  ABZeroPotential(std::shared_ptr<BasisController> basisA, std::shared_ptr<BasisController> basisB)
    : ABPotential<SCFMode>(basisA, basisB) {
  }
  /**
   * @brief Destructor
   */
  virtual ~ABZeroPotential() = default;

  /**
   * @brief Getter for the AB fock matrix contribution.
   * @return The AB fock matrix contribution.
   */
  SPMatrix<SCFMode>& getMatrix() {
    if (!_abPotential)
      _abPotential.reset(new SPMatrix<SCFMode>(this->_basisA->getNBasisFunctions(), this->_basisB->getNBasisFunctions()));
    return *_abPotential;
  }

 private:
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
};
template class ABZeroPotential<Options::SCF_MODES::RESTRICTED>;
template class ABZeroPotential<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABZEROPOTENTIAL_H_ */
