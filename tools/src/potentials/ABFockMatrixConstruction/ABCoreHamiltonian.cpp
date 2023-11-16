/**
 * @file ABCoreHamiltonian.cpp
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

/* Include Class Header*/
#include "potentials/ABFockMatrixConstruction/ABCoreHamiltonian.h"
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "potentials/ABFockMatrixConstruction/ABEffectiveCorePotential.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ABCoreHamiltonian<SCFMode>::ABCoreHamiltonian(std::shared_ptr<BasisController> basisA,
                                              std::shared_ptr<BasisController> basisB, std::shared_ptr<Geometry> geometry)
  : ABPotential<SCFMode>(basisA, basisB), _geom(geometry), _libint(Libint::getSharedPtr()) {
  // Setting the notifying system up.
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  _abEffectiveCorePotential =
      std::make_shared<ABEffectiveCorePotential<SCFMode>>(this->_basisA, this->_basisB, _geom->getAtoms());
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABCoreHamiltonian<SCFMode>::getMatrix() {
  if (!_abPotential) {
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    // initialize fock matrix
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA, nBasisB));
    SPMatrix<SCFMode>& f_AB = *_abPotential;
    Eigen::MatrixXd ints = _libint->compute1eInts(LIBINT_OPERATOR::kinetic, this->_basisA, this->_basisB).transpose();
    ints += _libint->compute1eInts(LIBINT_OPERATOR::nuclear, this->_basisA, this->_basisB, _geom->getAtoms()).transpose();
    for_spin(f_AB) {
      f_AB_spin = ints;
    };
    // add ECP contribution. (Zero if no ECPs are used)
    f_AB += _abEffectiveCorePotential->getMatrix();
  } /* if !_abPotential */
  return *_abPotential;
}

template class ABCoreHamiltonian<Options::SCF_MODES::RESTRICTED>;
template class ABCoreHamiltonian<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
