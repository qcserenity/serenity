/**
 * @file ElectronicStructureOptions.cpp
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
#include "settings/ElectronicStructureOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
template<>
void resolve<ORBITAL_FILE_TYPES>(std::string& value, ORBITAL_FILE_TYPES& field) {
  static const std::map<std::string, ORBITAL_FILE_TYPES> m = {{"SERENITY", ORBITAL_FILE_TYPES::SERENITY},
                                                              {"TURBOMOLE", ORBITAL_FILE_TYPES::TURBOMOLE},
                                                              {"MOLPRO", ORBITAL_FILE_TYPES::MOLPRO},
                                                              {"MOLCAS", ORBITAL_FILE_TYPES::MOLCAS},
                                                              {"MOLDEN", ORBITAL_FILE_TYPES::MOLDEN}};
  check(m, value, field);
}

template<>
void resolve<INTEGRAL_FILE_TYPES>(std::string& value, INTEGRAL_FILE_TYPES& field) {
  static const std::map<std::string, INTEGRAL_FILE_TYPES> m = {{"ASCII", INTEGRAL_FILE_TYPES::ASCII},
                                                               {"HDF5", INTEGRAL_FILE_TYPES::HDF5}};
  check(m, value, field);
}

template<>
void resolve<SCF_MODES>(std::string& value, SCF_MODES& field) {
  static const std::map<std::string, SCF_MODES> m = {{"RESTRICTED", SCF_MODES::RESTRICTED},
                                                     {"UNRESTRICTED", SCF_MODES::UNRESTRICTED}};
  check(m, value, field);
}

template<>
void resolve<ROHF_TYPES>(std::string& value, ROHF_TYPES& field) {
  static const std::map<std::string, ROHF_TYPES> m = {
      {"NONE", ROHF_TYPES::NONE}, {"CUHF", ROHF_TYPES::CUHF}, {"SUHF", ROHF_TYPES::SUHF}};
  check(m, value, field);
}

template<>
void resolve<ELECTRONIC_STRUCTURE_THEORIES>(std::string& value, ELECTRONIC_STRUCTURE_THEORIES& field) {
  static const std::map<std::string, ELECTRONIC_STRUCTURE_THEORIES> m = {{"HF", ELECTRONIC_STRUCTURE_THEORIES::HF},
                                                                         {"DFT", ELECTRONIC_STRUCTURE_THEORIES::DFT}};
  check(m, value, field);
}

template<>
void resolve<GLOBAL_PRINT_LEVELS>(std::string& value, GLOBAL_PRINT_LEVELS& field) {
  static const std::map<std::string, GLOBAL_PRINT_LEVELS> m = {
      {"MINIMUM", GLOBAL_PRINT_LEVELS::MINIMUM},     {"0", GLOBAL_PRINT_LEVELS::MINIMUM},
      {"NORMAL", GLOBAL_PRINT_LEVELS::NORMAL},       {"1", GLOBAL_PRINT_LEVELS::NORMAL},
      {"VERBOSE", GLOBAL_PRINT_LEVELS::VERBOSE},     {"2", GLOBAL_PRINT_LEVELS::VERBOSE},
      {"DEBUG", GLOBAL_PRINT_LEVELS::DEBUGGING},     {"3", GLOBAL_PRINT_LEVELS::DEBUGGING},
      {"DEBUGGING", GLOBAL_PRINT_LEVELS::DEBUGGING}, {"3", GLOBAL_PRINT_LEVELS::DEBUGGING}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
