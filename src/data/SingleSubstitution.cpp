/**
 * @file SingleSubstitution.cpp
 *
 * @date May 23, 2019
 * @author Moritz Bensberg
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
#include "data/SingleSubstitution.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h" //Definition of an orbital pair.

namespace Serenity {

SingleSubstitution::SingleSubstitution(std::shared_ptr<OrbitalPair> diagonalPair, double singlesPNOScaling)
  : i(diagonalPair->i), diagPair(diagonalPair), _pnoThreshold(diagonalPair->getPNOThreshold() * singlesPNOScaling){};

double SingleSubstitution::getPNOThreshold() {
  return _pnoThreshold;
}

} /* namespace Serenity */
