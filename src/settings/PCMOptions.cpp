/**
 * @file PCMOptions.cpp
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
#include "settings/PCMOptions.h"
/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<PCM_CAVITY_TYPES>(std::string& value, PCM_CAVITY_TYPES& field) {
  static const std::map<std::string, PCM_CAVITY_TYPES> m = {{"GEPOL_SAS", PCM_CAVITY_TYPES::GEPOL_SAS},
                                                            {"GEPOL_SES", PCM_CAVITY_TYPES::GEPOL_SES},
                                                            {"DELLEY", PCM_CAVITY_TYPES::DELLEY}};
  check(m, value, field);
}

template<>
void resolve<PCM_ATOMIC_RADII_TYPES>(std::string& value, PCM_ATOMIC_RADII_TYPES& field) {
  static const std::map<std::string, PCM_ATOMIC_RADII_TYPES> m = {{"BONDI", PCM_ATOMIC_RADII_TYPES::BONDI},
                                                                  {"UFF", PCM_ATOMIC_RADII_TYPES::UFF}};
  check(m, value, field);
}

template<>
void resolve<PCM_SOLVENTS>(std::string& value, PCM_SOLVENTS& field) {
  static const std::map<std::string, PCM_SOLVENTS> m = {{"H2O", PCM_SOLVENTS::WATER},
                                                        {"WATER", PCM_SOLVENTS::WATER},
                                                        {"PROPYLENECARBONATE", PCM_SOLVENTS::PROPYLENE_CARBONATE},
                                                        {"C4H6O3", PCM_SOLVENTS::PROPYLENE_CARBONATE},
                                                        {"DMSO", PCM_SOLVENTS::DMSO},
                                                        {"DIMETHYLSULFOXIDE", PCM_SOLVENTS::DMSO},
                                                        {"C2H6OS", PCM_SOLVENTS::DMSO},
                                                        {"CH3NO2", PCM_SOLVENTS::NITROMETHANE},
                                                        {"NITROMETHANE", PCM_SOLVENTS::NITROMETHANE},
                                                        {"CH3CN", PCM_SOLVENTS::ACETONITRILE},
                                                        {"ACETONITRILE", PCM_SOLVENTS::ACETONITRILE},
                                                        {"CH3OH", PCM_SOLVENTS::METHANOL},
                                                        {"MEOH", PCM_SOLVENTS::METHANOL},
                                                        {"METHANOL", PCM_SOLVENTS::METHANOL},
                                                        {"CH3CH2OH", PCM_SOLVENTS::ETHANOL},
                                                        {"ETOH", PCM_SOLVENTS::ETHANOL},
                                                        {"ETHANOL", PCM_SOLVENTS::ETHANOL},
                                                        {"C2H6CO", PCM_SOLVENTS::ACETONE},
                                                        {"ACETONE", PCM_SOLVENTS::ACETONE},
                                                        {"1,2-DICHLORETHANE", PCM_SOLVENTS::DICHLORETHANE},
                                                        {"DICHLORETHANE", PCM_SOLVENTS::DICHLORETHANE},
                                                        {"C2H4CL2", PCM_SOLVENTS::DICHLORETHANE},
                                                        {"METHYLENECHLORIDE", PCM_SOLVENTS::METHYLENECHLORIDE},
                                                        {"CH2CL2", PCM_SOLVENTS::METHYLENECHLORIDE},
                                                        {"THF", PCM_SOLVENTS::THF},
                                                        {"TETRAHYDROFURANE", PCM_SOLVENTS::THF},
                                                        {"C4H8O", PCM_SOLVENTS::THF},
                                                        {"C6H5NH2", PCM_SOLVENTS::ANILINE},
                                                        {"ANILINE", PCM_SOLVENTS::ANILINE},
                                                        {"C6H5CL", PCM_SOLVENTS::CHLOROBENZENE},
                                                        {"CHLOROBENZENE", PCM_SOLVENTS::CHLOROBENZENE},
                                                        {"CHCL3", PCM_SOLVENTS::CHLOROFORM},
                                                        {"CHLOROFORM", PCM_SOLVENTS::CHLOROFORM},
                                                        {"C6H5CH3", PCM_SOLVENTS::TOLUENE},
                                                        {"TOLUENE", PCM_SOLVENTS::TOLUENE},
                                                        {"DIOXANE", PCM_SOLVENTS::DIOXANE},
                                                        {"C4H8O2", PCM_SOLVENTS::DIOXANE},
                                                        {"1,4-DIOXANE", PCM_SOLVENTS::DIOXANE},
                                                        {"C6H6", PCM_SOLVENTS::BENZENE},
                                                        {"BENZENE", PCM_SOLVENTS::BENZENE},
                                                        {"CARBONTETRACHLORIDE", PCM_SOLVENTS::CARBON_TETRACHLORIDE},
                                                        {"CCL4", PCM_SOLVENTS::CARBON_TETRACHLORIDE},
                                                        {"C6H12", PCM_SOLVENTS::CYCLOHEXANE},
                                                        {"CYCLOHEXANE", PCM_SOLVENTS::CYCLOHEXANE},
                                                        {"C7H16", PCM_SOLVENTS::N_HEPTANE},
                                                        {"N-HEPTANE", PCM_SOLVENTS::N_HEPTANE},
                                                        {"EXPLICIT", PCM_SOLVENTS::EXPLICIT}};
  check(m, value, field);
}

template<>
void resolve<PCM_SOLVER_TYPES>(std::string& value, PCM_SOLVER_TYPES& field) {
  static const std::map<std::string, PCM_SOLVER_TYPES> m = {{"IEFPCM", PCM_SOLVER_TYPES::IEFPCM},
                                                            {"CPCM", PCM_SOLVER_TYPES::CPCM}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
