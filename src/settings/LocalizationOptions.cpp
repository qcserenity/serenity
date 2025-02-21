/**
 * @file LocalizationOptions.cpp
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
#include "settings/LocalizationOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
template<>
void resolve<ORBITAL_LOCALIZATION_ALGORITHMS>(std::string& value, ORBITAL_LOCALIZATION_ALGORITHMS& field) {
  static const std::map<std::string, ORBITAL_LOCALIZATION_ALGORITHMS> m = {
      {"PM", ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY},
      {"FB", ORBITAL_LOCALIZATION_ALGORITHMS::FOSTER_BOYS},
      {"IAO", ORBITAL_LOCALIZATION_ALGORITHMS::IAO},
      {"IBO", ORBITAL_LOCALIZATION_ALGORITHMS::IBO},
      {"ER", ORBITAL_LOCALIZATION_ALGORITHMS::EDMISTON_RUEDENBERG},
      {"NO", ORBITAL_LOCALIZATION_ALGORITHMS::NON_ORTHOGONAL},
      {"ALIGN", ORBITAL_LOCALIZATION_ALGORITHMS::ALIGN},
      {"NONE", ORBITAL_LOCALIZATION_ALGORITHMS::NONE}};
  check(m, value, field);
}

template<>
void resolve<POPULATION_ANALYSIS_ALGORITHMS>(std::string& value, POPULATION_ANALYSIS_ALGORITHMS& field) {
  static const std::map<std::string, POPULATION_ANALYSIS_ALGORITHMS> m = {
      {"MUL", POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN}, {"HIRSHFELD", POPULATION_ANALYSIS_ALGORITHMS::HIRSHFELD},
      {"IAO", POPULATION_ANALYSIS_ALGORITHMS::IAO},      {"IAOSHELL", POPULATION_ANALYSIS_ALGORITHMS::IAOShell},
      {"BECKE", POPULATION_ANALYSIS_ALGORITHMS::BECKE},  {"CM5", POPULATION_ANALYSIS_ALGORITHMS::CM5},
      {"CHELPG", POPULATION_ANALYSIS_ALGORITHMS::CHELPG}};
  check(m, value, field);
}

template<>
void resolve<DOS_SETTINGS>(std::string& value, DOS_SETTINGS& field) {
  static const std::map<std::string, DOS_SETTINGS> m = {
      {"LOOSE", DOS_SETTINGS::LOOSE},           {"NORMAL", DOS_SETTINGS::NORMAL},
      {"TIGHT", DOS_SETTINGS::TIGHT},           {"VERYTIGHT", DOS_SETTINGS::VERY_TIGHT},
      {"VERY_TIGHT", DOS_SETTINGS::VERY_TIGHT}, {"EXTREME", DOS_SETTINGS::EXTREME},
      {"SPREAD", DOS_SETTINGS::SPREAD}};
  check(m, value, field);
}
} /* namespace Options */
} /* namespace Serenity */
