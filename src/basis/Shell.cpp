/**
 * @file Shell.cpp
 *
 * @date Nov 7, 2016
 * @author Jan Unsleber
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
#include "basis/Shell.h"
/* Include Serenity Internal Headers */
#include "integrals/Normalization.h"
/* Include Std and External Headers */
#include <boost/container/small_vector.hpp>

namespace Serenity {

Shell::Shell(libint2::svector<double> exponents, libint2::svector<double> contractions, unsigned int angularMomentum,
             bool spherical, std::array<double, 3> coords, std::string element)
  : libint2::Shell(exponents, {{(int)angularMomentum, spherical, contractions}}, coords),
    _normFactors(
        new Eigen::VectorXd(!(this->contr[0].pure) ? N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l])),
    _contractions(contractions),
    _exponents(exponents),
    _element(element) {
  if (!this->contr[0].pure) {
    const double undoLibint = 1.0 / Normalization::finalNormalization(this->contr[0].l, 0, 0);
    for (int aX = this->contr[0].l, index = 0; aX >= 0; --aX) {
      for (int aY = this->contr[0].l - aX; aY >= 0; --aY) {
        int aZ = this->contr[0].l - aY - aX;
        (*_normFactors)[index] = undoLibint * Normalization::finalNormalization(aX, aY, aZ);
        ++index;
      }
    }
  }
}

Shell::Shell(libint2::svector<double> exponents, libint2::svector<double> exponents1, libint2::svector<double> contractions,
             unsigned int angularMomentum, bool spherical, std::array<double, 3> coords, std::string element)
  : libint2::Shell(std::move(exponents), {{(int)angularMomentum, spherical, contractions}}, coords),
    _normFactors(
        new Eigen::VectorXd(!(this->contr[0].pure) ? N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l])),
    _contractions(contractions),
    _exponents(exponents1),
    _element(element) {
  if (!this->contr[0].pure) {
    const double undoLibint = 1.0 / Normalization::finalNormalization(this->contr[0].l, 0, 0);
    for (int aX = this->contr[0].l, index = 0; aX >= 0; --aX) {
      for (int aY = this->contr[0].l - aX; aY >= 0; --aY) {
        int aZ = this->contr[0].l - aY - aX;
        (*_normFactors)[index] = undoLibint * Normalization::finalNormalization(aX, aY, aZ);
        ++index;
      }
    }
  }
}

Shell::Shell(const Shell& other) : libint2::Shell(other), NotifyingClass<Shell>() {
  _normFactors.reset(
      new Eigen::VectorXd(!(this->contr[0].pure) ? N_SHELL_CART[this->contr[0].l] : N_SHELL_SPH[this->contr[0].l]));
  if (!this->contr[0].pure) {
    const double undoLibint = 1.0 / Normalization::finalNormalization(this->contr[0].l, 0, 0);
    for (int aX = this->contr[0].l, index = 0; aX >= 0; --aX) {
      for (int aY = this->contr[0].l - aX; aY >= 0; --aY) {
        int aZ = this->contr[0].l - aY - aX;
        (*_normFactors)[index] = undoLibint * Normalization::finalNormalization(aX, aY, aZ);
        ++index;
      }
    }
  }
}

bool Shell::operator==(const Shell& other) const {
  // Check origin
  bool diffO =
      (std::fabs(this->O[0] - other.O[0]) + std::fabs(this->O[1] - other.O[1]) + std::fabs(this->O[2] - other.O[2])) > 5e-6;
  if (diffO)
    return false;
  // Check angular momentum
  if (getAngularMomentum() != other.getAngularMomentum())
    return false;
  // Check exponents and contractions.
  const auto otherExponents = other.getExponents();
  const auto otherContractions = other.getContractions();
  bool diffExpCont = _exponents.size() == otherExponents.size();
  if (diffExpCont) {
    assert(otherExponents.size() == otherContractions.size());
    for (unsigned int iPrim = 0; iPrim < _exponents.size(); ++iPrim) {
      double diffExp = std::fabs(_exponents[iPrim] - otherExponents[iPrim]);
      double diffCont = std::fabs(_contractions[iPrim] - otherContractions[iPrim]);
      if (diffExp > 1e-12 || diffCont > 1e-12) {
        return false;
      } // if diffExp > 1e-12 || diffCont > 1e-12
    }   // for unsigned int iPrim=0; iPrim < _exponents.size();++iPrim
  }     // if diffExpCont
  else {
    return false;
  }
  // Check for spherical or cartesian function.
  if (isSpherical() != other.isSpherical())
    return false;
  // If all checks were fine, these shells are the same.
  return true;
};

} /* namespace Serenity */
