/**
 * @file BasisOptions.h
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

#ifndef SETTINGS_BASISOPTIONS_H_
#define SETTINGS_BASISOPTIONS_H_

#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                           Basis                                                */
/**************************************************************************************************/
/**
 * For what a certain basis set should be used.\n
 * DEFAULT:          The basis set in which the electronic structure is calculated and expressed\n
 * AUX_COULOMB:      The auxiliary basis for a density fitting to evaluate Coulombic electron-electron
 *                   interactions
 * MINBAS            Default minimal basis set.
 * HUECKEL:          The basis set used for semiempirical calculations, most probably a minimal basis
 * IAO_LOCALIZATION: The basis set in which the Intrinsic Atomic Orbitals will be expressed
 * SCF_DENS_GUESS:   Basis set for the atom density guess.
 * AUX_CORREL:       The auxiliary basis for density fitting to evaluate electron-electron correlation
 *                   contributions
 */
enum class BASIS_PURPOSES {
  DEFAULT = 0,
  AUX_COULOMB = 1,
  MINBAS = 2,
  HUECKEL = 3,
  IAO_LOCALIZATION = 4,
  SCF_DENS_GUESS = 5,
  AUX_CORREL = 6
};
template<>
void resolve<BASIS_PURPOSES>(std::string& value, BASIS_PURPOSES& field);

/**************************************************************************************************/
/*                                     Density Fitting                                            */
/**************************************************************************************************/
/**
 * How a Hartree-Potential should be calculated.
 * (Old explanation:
 * In case of a split Coulomb and exchange part in the Fock update (e.g. in DFT no Exchange is
 * needed) these are the possibilities to evaluate the Coulomb part. In principle all complete
 * HARTREE_FOCK_POTENTIAL_CALCULATORS are possible, but also the famous RI approximation, which is a
 * way to reduce the scaling behaviour of DFT from O(n^4) to O(n^3).)
 */
enum class DENS_FITS { RI = 0, NONE = 1 };
template<>
void resolve<DENS_FITS>(std::string& value, DENS_FITS& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_BASISOPTIONS_H_ */
