/**
 * @file CorrelatedMethodsOptions.cpp
 *
 * @author Moritz Bensberg
 * @date May 12, 2020
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
#include "settings/CorrelatedMethodsOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<CC_LEVEL>(std::string& value, CC_LEVEL& field) {
  static const std::map<std::string, CC_LEVEL> m = {{"CCSD", CC_LEVEL::CCSD},
                                                    {"CCSD(T)", CC_LEVEL::CCSD_T},
                                                    {"DLPNO-CCSD", CC_LEVEL::DLPNO_CCSD},
                                                    {"DLPNO-CCSD(T0)", CC_LEVEL::DLPNO_CCSD_T0}};
  check(m, value, field);
}

template<>
void resolve<MP2_TYPES>(std::string& value, MP2_TYPES& field) {
  static const std::map<std::string, MP2_TYPES> m = {
      {"AO", MP2_TYPES::AO}, {"DF", MP2_TYPES::DF}, {"LOCAL", MP2_TYPES::LOCAL}, {"LT", MP2_TYPES::LT}};
  check(m, value, field);
}

template<>
void resolve<PNO_SETTINGS>(std::string& value, PNO_SETTINGS& field) {
  static const std::map<std::string, PNO_SETTINGS> m = {
      {"LOOSE", PNO_SETTINGS::LOOSE}, {"NORMAL", PNO_SETTINGS::NORMAL}, {"TIGHT", PNO_SETTINGS::TIGHT}};
  check(m, value, field);
}

template<>
void resolve<PNO_METHOD>(std::string& value, PNO_METHOD& field) {
  static const std::map<std::string, PNO_METHOD> m = {{"DLPNO-MP2", PNO_METHOD::DLPNO_MP2},
                                                      {"DLPNO_MP2", PNO_METHOD::DLPNO_MP2},
                                                      {"LMP2", PNO_METHOD::DLPNO_MP2},
                                                      {"DLPNO-CCSD", PNO_METHOD::DLPNO_CCSD},
                                                      {"DLPNO-CCSD(T0)", PNO_METHOD::DLPNO_CCSD_T0},
                                                      {"SC-MP2", PNO_METHOD::SC_MP2},
                                                      {"NONE", PNO_METHOD::NONE},
                                                      {"HF", PNO_METHOD::NONE}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
