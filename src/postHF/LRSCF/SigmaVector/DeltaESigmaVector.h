/**
 * @file DeltaESigmaVector.h
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

#ifndef DELTAESIGMAVECTOR_H_
#define DELTAESIGMAVECTOR_H_

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/SigmaVector/SigmaVector.h"

namespace Serenity {

template<Options::SCF_MODES T> class DeltaESigmaVector : public SigmaVector<T> {
public:
  DeltaESigmaVector(
      std::shared_ptr<SystemController> system,
      Eigen::MatrixXd& guessVector);
  virtual ~DeltaESigmaVector() = default;
protected:
  Eigen::MatrixXd calculateBlock(
      Eigen::MatrixXd& guess,
      std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens);
};

} /* namespace Serenity */

#endif /* DELTAESIGMAVECTOR_H_ */
