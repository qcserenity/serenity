/**
 * @file   StaticDamping.cpp
 *
 * @date   29. Dezember 2013, 14:13
 * @author Thomas Dresselhaus
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
#include "scf/damper/StaticDamping.h"

namespace Serenity {
using namespace std;

template<Options::SCF_MODES T>
StaticDamping<T>::StaticDamping(const double dampingFactor)
  : Damper<T>(), _oldMatrix(), _dampingFactor(dampingFactor), _initialized(false) {
}

template<Options::SCF_MODES T>
void StaticDamping<T>::damp(FockMatrix<T>& newMatrix) {
  if (_initialized) {
    for_spin(_oldMatrix, newMatrix) {
      assert(_oldMatrix_spin.cols() == newMatrix_spin.cols());
      assert(_oldMatrix_spin.rows() == newMatrix_spin.rows());
      newMatrix_spin *= (1.0 - _dampingFactor);
      newMatrix_spin += (_oldMatrix_spin * _dampingFactor);
    };
  }
  /*
   * Store as the reference for the next cycle.
   */
  for_spin(_oldMatrix, newMatrix) {
    _oldMatrix_spin = newMatrix_spin;
  };
  _initialized = true;
}

template<Options::SCF_MODES T>
void StaticDamping<T>::damp(SpinPolarizedData<T, Eigen::MatrixXd>& newMatrix) {
  if (_initialized) {
    for_spin(_oldMatrix, newMatrix) {
      assert(_oldMatrix_spin.cols() == newMatrix_spin.cols());
      assert(_oldMatrix_spin.rows() == newMatrix_spin.rows());
      newMatrix_spin *= (1.0 - _dampingFactor);
      newMatrix_spin += _oldMatrix_spin * _dampingFactor;
    };
  }
  /*
   * Store as the reference for the next cycle.
   */
  for_spin(_oldMatrix, newMatrix) {
    _oldMatrix_spin = newMatrix_spin;
  };
  _initialized = true;
}

template class StaticDamping<Options::SCF_MODES::RESTRICTED>;
template class StaticDamping<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
