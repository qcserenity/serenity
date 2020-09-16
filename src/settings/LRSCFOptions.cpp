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

#include "settings/LRSCFOptions.h"
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<RESPONSE_PROBLEM>(std::string& value, RESPONSE_PROBLEM& field) {
  static const std::map<std::string, RESPONSE_PROBLEM> m = {
      {"TDA", RESPONSE_PROBLEM::TDA},
      {"TDDFT", RESPONSE_PROBLEM::TDDFT},
      {"RPA", RESPONSE_PROBLEM::RPA},
  };
  check(m, value, field);
}

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
void resolve<MULTIPLICITY>(std::string& value, MULTIPLICITY& field) {
  static const std::map<std::string, MULTIPLICITY> m = {{"SINGLET", MULTIPLICITY::SINGLET}, {"TRIPLET", MULTIPLICITY::TRIPLET}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
