/**
 * @file SCFOptions.h
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

#ifndef SETTINGS_SCFOPTIONS_H_
#define SETTINGS_SCFOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                            SCF                                                 */
/**************************************************************************************************/
/**
 * How the initial guess for electronic structure calculations is determined\n
 * H_CORE:    hCore guess, i.e. ignoring electron-electron interactions,
 *            points to HCoreGuessCalculator\n
 * EHT:       Extended Hueckel Theory guess, points to ExtendedHueckel\n
 * ATOM_SCF:  The atomic densities are calculated with a proper SCF for each atom type.
 * SAP:       Superposition of atomic potentials.
 */
enum class INITIAL_GUESSES { H_CORE = 0, EHT = 1, ATOM_SCF = 2, ATOM_SCF_INPLACE = 3, SAP = 4 };
template<>
void resolve<INITIAL_GUESSES>(std::string& value, INITIAL_GUESSES& field);
/**
 * Which algorithm to use for the damping method.\n
 * NONE:    No damping will be applied.\n
 * STATIC:  Static damping using a constant factor.\n
 * SERIES:  The damping factor is changed from iteration to iteration.
 * DYNAMIC: Dynamic adjustment of the damping factor. Not implemented yet.
 */
enum class DAMPING_ALGORITHMS { NONE = 0, STATIC = 1, SERIES = 2, DYNAMIC = 3 };
template<>
void resolve<DAMPING_ALGORITHMS>(std::string& value, DAMPING_ALGORITHMS& field);

} /* namespace Options */
} /* namespace Serenity */
#endif /* SETTINGS_SCFOPTIONS_H_ */
