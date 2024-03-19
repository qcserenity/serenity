/**
 * @file   Damper.cpp
 *
 * @date   Feb 07, 2024
 * @author Lukas Paetow
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
#include "scf/damper/Damper.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Damper<SCFMode>::Damper() : _oldFock(), _oldDensity(), _initialized(false) {
}

template<Options::SCF_MODES SCFMode>
Damper<SCFMode>::Damper(const double dampingFactor) : _oldFock(), _dampingFactor(dampingFactor), _initialized(false) {
}

template<Options::SCF_MODES SCFMode>
Damper<SCFMode>::Damper(const double dStart, const double dStep, const double dEnd, int iStartUp)
  : _dampingFactor(dStart), _initialized(false), _dStart(dStart), _dStep(dStep), _dEnd(dEnd), _iStartUp(iStartUp + 1) {
  if (!(_iStartUp >= 0 and _dStart > 0.0 and _dEnd > 0.0 and _dStep > 0.0 and _dStart >= _dEnd - _dStep)) {
    throw SerenityError("Arithmetic series damping parameters are not set up correctly.");
  }
}

template<Options::SCF_MODES SCFMode>
void Damper<SCFMode>::staticDamp(FockMatrix<SCFMode>& newFock) {
  if (_initialized) {
    for_spin(_oldFock, newFock) {
      if (_oldFock_spin.rows() != newFock_spin.rows() || _oldFock_spin.cols() != newFock_spin.cols()) {
        throw SerenityError("Damping: Fock matrix dimensions of previous and current iteration do not match.");
      }
      newFock_spin *= (1.0 - _dampingFactor);
      newFock_spin += (_oldFock_spin * _dampingFactor);
    };
  }
  /*
   * Store as the reference for the next cycle.
   */
  for_spin(_oldFock, newFock) {
    _oldFock_spin = newFock_spin;
  };
  _initialized = true;
}

template<Options::SCF_MODES SCFMode>
void Damper<SCFMode>::dynamicDamp(FockMatrix<SCFMode>& newFock, DensityMatrix<SCFMode> newDensity) {
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

template<Options::SCF_MODES SCFMode>
void Damper<SCFMode>::arithmeticSeriesDamp(FockMatrix<SCFMode>& newFock) {
  if (_initialized) {
    // Update iStartUp.
    _iStartUp -= 1;

    // Update damping factor if startup is finished
    if (_iStartUp <= 0) {
      _dampingFactor -= _dStep;
    }

    // If damping factor is smaller than desired end value, set damping factor to end value
    if (_dampingFactor < _dEnd)
      _dampingFactor = _dEnd;

    // Damp
    for_spin(_oldFock, newFock) {
      if (_oldFock_spin.rows() != newFock_spin.rows() || _oldFock_spin.cols() != newFock_spin.cols()) {
        throw SerenityError("Damping: Fock matrix dimensions of previous and current iteration do not match.");
      }
      newFock_spin *= (1.0 - _dampingFactor);
      newFock_spin += _oldFock_spin * _dampingFactor;
    };
  }

  // Store current matrix as reference for the next cycle.
  for_spin(_oldFock, newFock) {
    _oldFock_spin = newFock_spin;
  };

  _initialized = true;
}

template class Damper<Options::SCF_MODES::RESTRICTED>;
template class Damper<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
