/**
 * @file NTOCalculator.h
 *
 * @date Jul 13, 2017
 * @author L. Hellmann and J. Gießbach
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

#ifndef NTOCALCULATOR_H_
#define NTOCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "system/SystemController.h"


namespace Serenity {

template<Options::SCF_MODES T>
class NTOCalculator {
public:

  NTOCalculator(
      std::shared_ptr<SystemController> activeSystem,
      const double plottingThreshold);

  virtual ~NTOCalculator() = default;

  /**
   * Returns the number of states stored by the LRSCF
   */
  unsigned int getNumberOfStates();

  std::string getDir(unsigned int i);

  SpinPolarizedData<T,Eigen::MatrixXd>& getOccNTOs(unsigned int iState) {
    if (!_hasBeenCalculated) {
      calcNTOs(iState);
    } else if (_state != iState) {
      calcNTOs(iState); 
    }
    return _occNTOs; 
  }
  SpinPolarizedData<T,Eigen::MatrixXd>& getVirtNTOs(unsigned int iState) {
    if (!_hasBeenCalculated) { 
      calcNTOs(iState);
    } else if (_state != iState) {
      calcNTOs(iState);
    }
    return _virtNTOs;
  }
  SpinPolarizedData<T,Eigen::VectorXd>& getOccEigenvalues(unsigned int iState) {
    if (!_hasBeenCalculated) { 
      calcNTOs(iState);
    } else if (_state != iState) {
      calcNTOs(iState);
    }
    return _occEigenvalues;
  }
  SpinPolarizedData<T,Eigen::VectorXd>& getVirtEigenvalues(unsigned int iState) {
    if (!_hasBeenCalculated) { 
      calcNTOs(iState);
    } else if (_state != iState) {
      calcNTOs(iState);
    }
    return _virtEigenvalues;
  }




private:

  void printNTOInfo(
      std::string spin,
      const unsigned int iState,
      Eigen::VectorXd& oEigenvalues,
      Eigen::VectorXd& vEigenvalues,
      Eigen::MatrixXd& u,
      Eigen::MatrixXd& v);

  bool readFromHDF5(
      Eigen::MatrixXd& X,
      Eigen::MatrixXd& Y,
      Eigen::VectorXd& eigenvalues);

  std::shared_ptr<SystemController> _activeSystem;
  const double _plottingThreshold;

  Eigen::MatrixXd _XPY;
  Eigen::VectorXd _eigenvalues;


  SpinPolarizedData<T, Eigen::MatrixXd> _occNTOs;
  SpinPolarizedData<T, Eigen::MatrixXd> _virtNTOs;
  SpinPolarizedData<T, Eigen::VectorXd> _occEigenvalues;
  SpinPolarizedData<T, Eigen::VectorXd> _virtEigenvalues;

  void calcNTOs(int iState);

  bool _hasBeenCalculated;
  unsigned int _state;

};

} /* namespace Serenity */

#endif /* NTOCALCULATOR_H_ */
