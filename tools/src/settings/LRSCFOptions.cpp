/**
 * @file LRSCFOptions.cpp
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
#include "settings/LRSCFOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<LRSCF_TYPE>(std::string& value, LRSCF_TYPE& field) {
  static const std::map<std::string, LRSCF_TYPE> m = {
      {"ISOLATED", LRSCF_TYPE::ISOLATED}, {"ISO", LRSCF_TYPE::ISOLATED},    {"UNCOUPLED", LRSCF_TYPE::UNCOUPLED},
      {"FDEU", LRSCF_TYPE::UNCOUPLED},    {"COUPLED", LRSCF_TYPE::COUPLED}, {"FDEC", LRSCF_TYPE::COUPLED},
  };
  check(m, value, field);
}

template<>
void resolve<INTEGRAL_TYPE>(std::string& value, INTEGRAL_TYPE& field) {
  static const std::map<std::string, INTEGRAL_TYPE> m = {{"NUMERICAL", INTEGRAL_TYPE::NUMERICAL},
                                                         {"ANALYTICAL", INTEGRAL_TYPE::ANALYTICAL}};
  check(m, value, field);
}

template<>
void resolve<GAUGE>(std::string& value, GAUGE& field) {
  static const std::map<std::string, GAUGE> m = {{"LENGTH", GAUGE::LENGTH}, {"VELOCITY", GAUGE::VELOCITY}};
  check(m, value, field);
}

template<>
void resolve<LR_METHOD>(std::string& value, LR_METHOD& field) {
  static const std::map<std::string, LR_METHOD> m = {
      {"CIS", LR_METHOD::TDA}, {"TDA", LR_METHOD::TDA},         {"RPA", LR_METHOD::TDDFT}, {"TDDFT", LR_METHOD::TDDFT},
      {"CC2", LR_METHOD::CC2}, {"CISDINF", LR_METHOD::CISDINF}, {"CISD", LR_METHOD::CISD}, {"ADC2", LR_METHOD::ADC2}};
  check(m, value, field);
}

template<>
void resolve<STABILITY_ANALYSIS>(std::string& value, STABILITY_ANALYSIS& field) {
  static const std::map<std::string, STABILITY_ANALYSIS> m = {{"NONE", STABILITY_ANALYSIS::NONE},
                                                              {"REAL", STABILITY_ANALYSIS::REAL},
                                                              {"NONREAL", STABILITY_ANALYSIS::NONREAL},
                                                              {"SPINFLIP", STABILITY_ANALYSIS::SPINFLIP}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
