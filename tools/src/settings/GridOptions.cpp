/**
 * @file GridOptions.cpp
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
#include "settings/GridOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<RADIAL_GRID_TYPES>(std::string& value, RADIAL_GRID_TYPES& field) {
  static const std::map<std::string, RADIAL_GRID_TYPES> m = {{"BECKE", RADIAL_GRID_TYPES::BECKE},
                                                             {"HANDY", RADIAL_GRID_TYPES::HANDY},
                                                             {"AHLRICHS", RADIAL_GRID_TYPES::AHLRICHS},
                                                             {"KNOWLES", RADIAL_GRID_TYPES::KNOWLES},
                                                             {"EQUIDISTAND", RADIAL_GRID_TYPES::EQUI}};
  check(m, value, field);
}

template<>
void resolve<SPHERICAL_GRID_TYPES>(std::string& value, SPHERICAL_GRID_TYPES& field) {
  static const std::map<std::string, SPHERICAL_GRID_TYPES> m = {{"LEBEDEV", SPHERICAL_GRID_TYPES::LEBEDEV}};
  check(m, value, field);
}

template<>
void resolve<GRID_TYPES>(std::string& value, GRID_TYPES& field) {
  static const std::map<std::string, GRID_TYPES> m = {
      {"BECKE", GRID_TYPES::BECKE}, {"VORONOI", GRID_TYPES::VORONOI}, {"SSF", GRID_TYPES::SSF}};
  check(m, value, field);
}

template<>
void resolve<GRID_PURPOSES>(std::string& value, GRID_PURPOSES& field) {
  static const std::map<std::string, GRID_PURPOSES> m = {
      {"DEFAULT", GRID_PURPOSES::DEFAULT}, {"SMALL", GRID_PURPOSES::SMALL}, {"PLOT", GRID_PURPOSES::PLOT}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
