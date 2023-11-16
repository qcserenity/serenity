/**
 * @file OrthogonalizationOptions.cpp
 *
 * @author Anja Massolle
 * @date Aug 06, 2020
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
#include "settings/OrthogonalizationOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
template<>
void resolve<ORTHOGONALIZATION_ALGORITHMS>(std::string& value, ORTHOGONALIZATION_ALGORITHMS& field) {
  static const std::map<std::string, ORTHOGONALIZATION_ALGORITHMS> m = {{"NONE", ORTHOGONALIZATION_ALGORITHMS::NONE},
                                                                        {"LOEWDIN", ORTHOGONALIZATION_ALGORITHMS::LOEWDIN},
                                                                        {"PIPEK", ORTHOGONALIZATION_ALGORITHMS::PIPEK},
                                                                        {"BROER", ORTHOGONALIZATION_ALGORITHMS::BROER}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */