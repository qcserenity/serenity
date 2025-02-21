/**
 * @file MBPTOptions.cpp
 *
 * @author Johannes Toelle
 * @date May 27, 2020
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
#include "settings/MBPTOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<MBPT>(std::string& value, MBPT& field) {
  static const std::map<std::string, MBPT> m = {
      {"GW", MBPT::GW},
      {"RPA", MBPT::RPA},
  };
  check(m, value, field);
}

template<>
void resolve<GWALGORITHM>(std::string& value, GWALGORITHM& field) {
  static const std::map<std::string, GWALGORITHM> m = {
      {"CD", GWALGORITHM::CD},
      {"AC", GWALGORITHM::AC},
      {"ANALYTIC", GWALGORITHM::ANALYTIC},
  };
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */