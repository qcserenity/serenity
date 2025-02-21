/**
 * @file GeometryOptions.cpp
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
#include "settings/GeometryOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<OPTIMIZATION_ALGORITHMS>(std::string& value, OPTIMIZATION_ALGORITHMS& field) {
  static const std::map<std::string, OPTIMIZATION_ALGORITHMS> m = {{"SQNM", OPTIMIZATION_ALGORITHMS::SQNM},
                                                                   {"BFGS", OPTIMIZATION_ALGORITHMS::BFGS}};
  check(m, value, field);
}

template<>
void resolve<GRADIENT_TYPES>(std::string& value, GRADIENT_TYPES& field) {
  static const std::map<std::string, GRADIENT_TYPES> m = {{"NUMERICAL", GRADIENT_TYPES::NUMERICAL},
                                                          {"ANALYTICAL", GRADIENT_TYPES::ANALYTICAL}};
  check(m, value, field);
}

template<>
void resolve<GEOMETRY_OPTIMIZATION_TYPES>(std::string& value, GEOMETRY_OPTIMIZATION_TYPES& field) {
  static const std::map<std::string, GEOMETRY_OPTIMIZATION_TYPES> m = {
      {"GROUNDSTATE", GEOMETRY_OPTIMIZATION_TYPES::GROUNDSTATE}, {"TS", GEOMETRY_OPTIMIZATION_TYPES::TS}};
  check(m, value, field);
}

template<>
void resolve<HESSIAN_TYPES>(std::string& value, HESSIAN_TYPES& field) {
  static const std::map<std::string, HESSIAN_TYPES> m = {{"NUMERICAL", HESSIAN_TYPES::NUMERICAL},
                                                         {"ANALYTICAL", HESSIAN_TYPES::ANALYTICAL}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
