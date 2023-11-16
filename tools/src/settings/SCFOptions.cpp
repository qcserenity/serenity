/**
 * @file SCFOptions.cpp
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
#include "settings/SCFOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<INITIAL_GUESSES>(std::string& value, INITIAL_GUESSES& field) {
  static const std::map<std::string, INITIAL_GUESSES> m = {{"HCORE", INITIAL_GUESSES::H_CORE},
                                                           {"EHT", INITIAL_GUESSES::EHT},
                                                           {"ATOM_SCF", INITIAL_GUESSES::ATOM_SCF},
                                                           {"ATOMSCF", INITIAL_GUESSES::ATOM_SCF},
                                                           {"ATOM_SCF_INPLACE", INITIAL_GUESSES::ATOM_SCF_INPLACE},
                                                           {"ATOMSCF_INPLACE", INITIAL_GUESSES::ATOM_SCF_INPLACE},
                                                           {"INPLACE", INITIAL_GUESSES::ATOM_SCF_INPLACE},
                                                           {"SAP", INITIAL_GUESSES::SAP}};
  check(m, value, field);
}

template<>
void resolve<DAMPING_ALGORITHMS>(std::string& value, DAMPING_ALGORITHMS& field) {
  static const std::map<std::string, DAMPING_ALGORITHMS> m = {{"NONE", DAMPING_ALGORITHMS::NONE},
                                                              {"STATIC", DAMPING_ALGORITHMS::STATIC},
                                                              {"SERIES", DAMPING_ALGORITHMS::SERIES},
                                                              {"DYNAMIC", DAMPING_ALGORITHMS::DYNAMIC}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
