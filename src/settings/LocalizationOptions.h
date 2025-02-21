/**
 * @file LocalizationOptions.h
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

#ifndef SETTINGS_LOCALIZATIONOPTIONS_H_
#define SETTINGS_LOCALIZATIONOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                       Localization                                             */
/**************************************************************************************************/

/**
 * Orbital localization methods.
 *  IAO:                 Project onto Intrinsic atomic orbital basis according to J. Chem. Theory Comput. 9, 483 (2013).
 *  IBO:                 Intrinsic bond orbital localization according to J. Chem. Theory Comput. 9, 4834 (2013).
 *  PIPEK_MEZEY:         Pipek--Mezey localization according to J. Chem. Phys. 90, 4916 (1989)
 *  FOSTER_BOYS:         Foster--Boys type localization according to J. Chem. Phys. 61, 3905 (1974)
 *  EDMISTON_RUEDENBERG: According to Rev. Mod. Phys. 35, 457 (1963)
 *  NON_ORTHOGONAL:      According to J. Chem. Phys. 120, 9458 (2004). This can lead to non-orthogonal orbitals.
 *  NONE:                No localization criterion.
 */
enum class ORBITAL_LOCALIZATION_ALGORITHMS {
  PIPEK_MEZEY = 0,
  FOSTER_BOYS = 1,
  IAO = 2,
  IBO = 3,
  EDMISTON_RUEDENBERG = 4,
  NON_ORTHOGONAL = 5,
  ALIGN = 6,
  NONE = 7
};
template<>
void resolve<ORBITAL_LOCALIZATION_ALGORITHMS>(std::string& value, ORBITAL_LOCALIZATION_ALGORITHMS& field);

/**************************************************************************************************/
/*                                   Population Analysis                                          */
/**************************************************************************************************/
enum class POPULATION_ANALYSIS_ALGORITHMS {
  MULLIKEN = 0,
  HIRSHFELD = 1,
  IAO = 2,
  IAOShell = 3,
  BECKE = 4,
  CM5 = 5,
  CHELPG = 6
};
template<>
void resolve<POPULATION_ANALYSIS_ALGORITHMS>(std::string& value, POPULATION_ANALYSIS_ALGORITHMS& field);

/**************************************************************************************************/
/*                                     DOS Macro Flags                                            */
/**************************************************************************************************/
/*
 * These flags can be used to automatically select a set of DOS thresholds. The accuracy of the
 * DOS orbital sets for relative energies is indicated by the flags name.
 */
enum class DOS_SETTINGS { LOOSE = 0, NORMAL = 1, TIGHT = 2, VERY_TIGHT = 3, EXTREME = 4, SPREAD = 5 };
template<>
void resolve<DOS_SETTINGS>(std::string& value, DOS_SETTINGS& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_LOCALIZATIONOPTIONS_H_ */
