/**
 * @file ABEffectiveCorePotential.cpp
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

/* Include Class Header*/
#include "potentials/ABFockMatrixConstruction/ABEffectiveCorePotential.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "integrals/wrappers/Libecpint.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABEffectiveCorePotential<SCFMode>::ABEffectiveCorePotential(std::shared_ptr<BasisController> basisA,
                                                            std::shared_ptr<BasisController> basisB,
                                                            std::vector<std::shared_ptr<Atom>> atoms)
  : ABPotential<SCFMode>(basisA, basisB), _atoms(atoms) {
  // Setting the notifying system up.
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  // Checking whether ECPs are actually present.
  for (const auto& atom : _atoms) {
    if (atom->usesECP())
      _notZero = true;
  }
}
template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABEffectiveCorePotential<SCFMode>::getMatrix() {
  if (!_abPotential) {
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    // initialize fock matrix
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA, nBasisB));
    if (_notZero) {
      // calculate the ECP integrals
      const auto ecpIntegrals = Libecpint::computeECPIntegrals(this->_basisA, this->_basisB, _atoms);
      // assign integrals
      auto& abPot = *_abPotential;
      for_spin(abPot) {
        abPot_spin = ecpIntegrals;
      };
    } // if _notZero
  }   // if !_abPotential
  return *_abPotential;
}

template class ABEffectiveCorePotential<Options::SCF_MODES::RESTRICTED>;
template class ABEffectiveCorePotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
