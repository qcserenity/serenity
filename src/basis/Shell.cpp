/**
 * @file Shell.cpp
 *
 * @date Nov 7, 2016
 * @author Jan Unsleber
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
#include "basis/Shell.h"
/* Include Serenity Internal Headers */
#include "integrals/Normalization.h"

namespace Serenity {

Shell::Shell(std::vector<double> exponents,
    std::vector<double> contractions,
    unsigned int angularMomentum,
    bool spherical,
    std::array<double,3> coords,
    std::string element ):
              libint2::Shell(exponents,
                  {{(int)angularMomentum,spherical,contractions}},
                  coords),
                  _normFactors(new Eigen::VectorXd(!(this->contr[0].pure) ?
                      N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l])),
                      _contractions(contractions),
                      _exponents(exponents),
                      _element(element){

  if (!this->contr[0].pure){
    const double undoLibint = 1.0/ Normalization::finalNormalization(this->contr[0].l, 0, 0);
    for (int aX = this->contr[0].l, index = 0; aX >= 0; --aX) {
      for (int aY = this->contr[0].l - aX; aY >= 0; --aY) {
        int aZ = this->contr[0].l - aY - aX;
        (*_normFactors)[index] = undoLibint * Normalization::finalNormalization(aX, aY, aZ);
        ++index;
      }
    }
  }
}

Shell::Shell(const Shell& other): libint2::Shell(other),
                                  NotifyingClass<Shell>(){
  _normFactors.reset(new Eigen::VectorXd(!(this->contr[0].pure) ?
      N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l]));
  if (!this->contr[0].pure){
    const double undoLibint = 1.0/ Normalization::finalNormalization(this->contr[0].l, 0, 0);
    for (int aX = this->contr[0].l, index = 0; aX >= 0; --aX) {
      for (int aY = this->contr[0].l - aX; aY >= 0; --aY) {
        int aZ = this->contr[0].l - aY - aX;
        (*_normFactors)[index] = undoLibint * Normalization::finalNormalization(aX, aY, aZ);
        ++index;
      }
    }
  }
}






} /* namespace Serenity */
