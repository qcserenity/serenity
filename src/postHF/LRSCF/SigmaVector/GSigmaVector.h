/**
 * @file GSigmaVector.h
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

#ifndef GSIGMAVECTOR_H_
#define GSIGMAVECTOR_H_

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/SigmaVector/SigmaVector.h"

namespace Serenity {

template<Options::SCF_MODES T>
class GSigmaVector: public SigmaVector<T> {
public:
  /**
   * @brief Sigma Vector from the (A+B) matrix in TDHF for singlet excitations.
   *        Calculates the second term of equation (11) in J. Chem. Phys 99, 1262, 1993.
   *        Note that there is a sign error in equation (14) defining G.
   *        In the case of TDDFT, exchange and long-range exchange are scaled with a
   *        scaling factor defined by the exchange-correlation functional.
   *
   * @param system
   * @param guessVector
   */
  GSigmaVector(
      std::shared_ptr<SystemController> system,
      Eigen::MatrixXd& guessVector);
  virtual ~GSigmaVector() = default;

protected:
  /**
   * @brief Calculates a block of sigma vectors.
   * @param guess The guess vectors.
   * @param dens The pseudo densities for each guess vector.
   * @return Returns the sigma vector block.
   */
  Eigen::MatrixXd calculateBlock(
      Eigen::MatrixXd& guess,
      std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens);

private:
  /// @brief The HF exchange ratio.
  double _hfExchangeRatio;
  /// @brief The long-range HF exchange ratio.
  double _lrhfExchangeRatio;
  /// @brief The range separation parameter.
  double _mu;
};

} /* namespace Serenity */

#endif /* GSIGMAVECTOR_H */
