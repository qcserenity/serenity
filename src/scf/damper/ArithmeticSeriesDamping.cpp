/**
 * @file ArithmeticSeriesDamping.cpp
 *
 * @date Oct 28, 2016
 * @author M. Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "scf/damper/ArithmeticSeriesDamping.h"



namespace Serenity {
template<Options::SCF_MODES T> ArithmeticSeriesDamping<T>::ArithmeticSeriesDamping(
    const double dStart,
    const double dStep,
    const double dEnd,
    int iStartUp):
     Damper<T>(),
     _dStart(dStart),
     _dStep(dStep),
     _dEnd(dEnd),
     _iStartUp(iStartUp + 1),
     _dampingFactor(dStart),
     _initialized(false){
  assert (_iStartUp >= 0 and _dStart > 0.0 and _dEnd > 0.0 and _dStep > 0.0);
  assert (_dStart >= _dEnd - _dStep);
}

template<Options::SCF_MODES T> void ArithmeticSeriesDamping<T>::damp(FockMatrix<T>& newMatrix) {
  if (_initialized){
    //Update iStartUp.
    _iStartUp -= 1;

    //Update damping factor if startup is finished
    if (_iStartUp <= 0) {
      _dampingFactor -= _dStep;
    }

    //If damping factor is smaller than desired end value, set damping factor to end value
    if (_dampingFactor < _dEnd ) _dampingFactor = _dEnd;

    //Damp
    for_spin(_oldMatrix,newMatrix){
      assert(_oldMatrix_spin.cols()==newMatrix_spin.cols());
      assert(_oldMatrix_spin.rows()==newMatrix_spin.rows());
      newMatrix_spin *= (1.0 - _dampingFactor);
      newMatrix_spin += _oldMatrix_spin * _dampingFactor;
    };
  }

  //Store current matrix as reference for the next cycle.
  for_spin(_oldMatrix,newMatrix){
    _oldMatrix_spin = newMatrix_spin;
  };

  _initialized = true;
}

template<Options::SCF_MODES T> void ArithmeticSeriesDamping<T>::damp(SpinPolarizedData<T, Eigen::MatrixXd>& newMatrix) {
  if (_initialized){
    //Update iStartUp.
    _iStartUp -= 1;

    //Update damping factor if startup is finished
    if (_iStartUp <= 0) {
      _dampingFactor -= _dStep;
    }

    //If damping factor is smaller than desired end value, set damping factor to end value
    if (_dampingFactor < _dEnd ) _dampingFactor = _dEnd;

    //Damp
    for_spin(_oldMatrix,newMatrix){
      assert(_oldMatrix_spin.cols()==newMatrix_spin.cols());
      assert(_oldMatrix_spin.rows()==newMatrix_spin.rows());
      newMatrix_spin *= (1.0 - _dampingFactor);
      newMatrix_spin += _oldMatrix_spin * _dampingFactor;
    };
  }

  //Store current matrix as reference for the next cycle.
  for_spin(_oldMatrix,newMatrix){
    _oldMatrix_spin = newMatrix_spin;
  };

  _initialized = true;
}


template class ArithmeticSeriesDamping<Options::SCF_MODES::RESTRICTED>;
template class ArithmeticSeriesDamping<Options::SCF_MODES::UNRESTRICTED>;

}/* namespace Serenity */

