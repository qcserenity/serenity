/**
 * @file PCMOptions.h
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

#ifndef SETTINGS_PCMOPTIONS_H_
#define SETTINGS_PCMOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

/**
 * Cavity types.
 */
enum class PCM_CAVITY_TYPES { GEPOL_SAS, GEPOL_SES, DELLEY };
template<>
void resolve<PCM_CAVITY_TYPES>(std::string& value, PCM_CAVITY_TYPES& field);

/**
 * Atomic radii types used in PCMSolver.
 */
enum class PCM_ATOMIC_RADII_TYPES { BONDI, UFF };
template<>
void resolve<PCM_ATOMIC_RADII_TYPES>(std::string& value, PCM_ATOMIC_RADII_TYPES& field);

/**
 * PCM solvents available. Note that not all solvents are available
 * since the interface to PCMSolver accepts solvent names up to a length
 * of 16 chars. Thus, all solvents with longer names have to be addressed
 * by selecting the "EXPLICIT" option and assigning all values manually.
 */
enum class PCM_SOLVENTS {
  WATER,
  PROPYLENE_CARBONATE,
  DMSO,
  NITROMETHANE,
  ACETONITRILE,
  METHANOL,
  ETHANOL,
  ACETONE,
  DICHLORETHANE,
  METHYLENECHLORIDE,
  THF,
  ANILINE,
  CHLOROBENZENE,
  CHLOROFORM,
  TOLUENE,
  DIOXANE,
  BENZENE,
  CARBON_TETRACHLORIDE,
  CYCLOHEXANE,
  N_HEPTANE,
  EXPLICIT
};
template<>
void resolve<PCM_SOLVENTS>(std::string& value, PCM_SOLVENTS& field);

/**
 * Solver used by PCMSolver.
 */
enum class PCM_SOLVER_TYPES { IEFPCM, CPCM };
template<>
void resolve<PCM_SOLVER_TYPES>(std::string& value, PCM_SOLVER_TYPES& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_PCMOPTIONS_H_ */
