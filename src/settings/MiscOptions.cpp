/**
 * @file MiscOptions.cpp
 *
 * @author Moritz Bensberg
 * @date May 11, 2020
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
#include "settings/MiscOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<SYSTEM_SPLITTING_ALGORITHM>(std::string& value, SYSTEM_SPLITTING_ALGORITHM& field) {
  static const std::map<std::string, SYSTEM_SPLITTING_ALGORITHM> m = {
      {"ENFORCE_CHARGES", SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES},
      {"ENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES},
      {"BEST_MATCH", SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH},
      {"BESTMATCH", SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH},
      {"POPULATION_THRESHOLD", SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD},
      {"POPULATION", SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD},
      {"SPADE", SYSTEM_SPLITTING_ALGORITHM::SPADE},
      {"SPADE_ENFORCE_CHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES},
      {"SPADE_ENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES},
      {"SPADEENFORCECHARGES", SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES}};
  check(m, value, field);
}

template<>
void resolve<BASIS_SET_TRUNCATION_ALGORITHMS>(std::string& value, BASIS_SET_TRUNCATION_ALGORITHMS& field) {
  static const std::map<std::string, BASIS_SET_TRUNCATION_ALGORITHMS> m = {
      {"NONE", BASIS_SET_TRUNCATION_ALGORITHMS::NONE},
      {"NETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION},
      {"NETPOPULATION", BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION},
      {"PRIMITIVENETPOP", BASIS_SET_TRUNCATION_ALGORITHMS::PRIMITIVE_NET_POPULATION},
  };
  check(m, value, field);
}

template<>
void resolve<GAUGE_ORIGIN>(std::string& value, GAUGE_ORIGIN& field) {
  static const std::map<std::string, GAUGE_ORIGIN> m = {{"CENTEROFMASS", GAUGE_ORIGIN::COM}, {"ORIGIN", GAUGE_ORIGIN::ORIGIN}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
