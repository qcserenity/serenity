/**
 * @file SigmaVector.h
 *
 * @date Oct 09, 2017
 * @author Michael Boeckers
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

#ifndef SIGMAVECTOR_H_
#define SIGMAVECTOR_H_


/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "data/SpinPolarizedData.h"
#include "system/SystemController.h"
#include "data/matrices/SPMatrix.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES T> class SigmaVector {
public:
  SigmaVector();
  SigmaVector(
      std::shared_ptr<SystemController> system,
      Eigen::MatrixXd& guessVector);
  virtual ~SigmaVector() = default;

  Eigen::MatrixXd& getSigma() {
    if (!_hasBeenCalculated) calculate();
    return _sigmaVector;
  }

protected:
  /**
   * @brief Calculates a block of sigma vectors.
   * @param guess The guess vectors.
   * @param dens The pseudo densities for each subsystem.
   * @return Returns the sigma vector block.
   */
  virtual Eigen::MatrixXd calculateBlock(
      Eigen::MatrixXd& guess,
      std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens) = 0;
  Eigen::VectorXd ao2mo(SPMatrix<T>& pF);
  SpinPolarizedData<T,Eigen::MatrixXd> pDens(unsigned int iGuess);
  std::shared_ptr<SystemController> _system;
  Eigen::MatrixXd _guessVector;
  Eigen::MatrixXd _sigmaVector;
  bool _hasBeenCalculated;
  void calculate();


};

} /* namespace Serenity */

#endif /* SIGMAVECTOR_H_ */
