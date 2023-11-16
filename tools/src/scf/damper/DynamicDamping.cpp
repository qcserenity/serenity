/**
 * @file   DynamicDamping.cpp
 *
 * @date   13 August 2020
 * @author Niklas Niemeyer
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
#include "scf/damper/DynamicDamping.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DynamicDamping<SCFMode>::DynamicDamping() : Damper<SCFMode>(), _oldFock(), _oldDensity(), _initialized(false) {
}

template<Options::SCF_MODES SCFMode>
void DynamicDamping<SCFMode>::dynamicDamp(FockMatrix<SCFMode>& newFock, DensityMatrix<SCFMode> newDensity) {
  if (_initialized) {
    for_spin(_oldFock, _oldDensity, newFock, newDensity) {
      double s = _oldFock_spin.cwiseProduct(newDensity_spin - _oldDensity_spin).sum();
      double c = (newFock_spin - _oldFock_spin).cwiseProduct(newDensity_spin - _oldDensity_spin).sum();
      double dampingFactor = (c <= -s / 2) ? 1.0 : -s / (2 * c);
      printf("   *** Dynamic damping factor: %10.3f ***\n", 1 - dampingFactor);

      newFock_spin = ((1 - dampingFactor) * _oldFock_spin + dampingFactor * newFock_spin).eval();
      newDensity_spin = ((1 - dampingFactor) * _oldDensity_spin + dampingFactor * newDensity_spin).eval();
    };
  }

  // store updated Fock and density matrices
  for_spin(_oldFock, _oldDensity, newFock, newDensity) {
    _oldFock_spin = newFock_spin;
    _oldDensity_spin = newDensity_spin;
  };

  _initialized = true;
}

template class DynamicDamping<Options::SCF_MODES::RESTRICTED>;
template class DynamicDamping<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
