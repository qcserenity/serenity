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

/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <string>

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
 * HUECKEL:          The basis set used for semiempirical calculations, most probably a minimal basis.
 * IAO_LOCALIZATION: The basis set in which the Intrinsic Atomic Orbitals will be expressed.
 * SCF_DENS_GUESS:   Basis set for the atom density guess.
 * AUX_CORREL:       The auxiliary basis for density fitting to evaluate electron-electron correlation
 *                   contributions.
 * ATOMIC_CHOLESKY:  The universal auxiliary basis generated using the atomic Cholesky decomposition.
 * ATOMIC_COMPACT_CHOLESKY:  The universal auxiliary basis generated using the atomic-compact Cholesky decomposition.
 * ERF_ATOMIC_CHOLESKY: The auxiliary basis generated for long-range contributions using the atomic Cholesky
 *                    decomposition.
 * ERF_ATOMIC_COMPACT_CHOLESKY: The auxiliary basis generated for long-range contributions using the
 *                    atomic-compact Cholesky decomposition.
 * AUX_JK:  The auxiliary basis for density fitting to evaluate Coulomb and exchange contributions.
 */
enum class BASIS_PURPOSES {
  DEFAULT = 0,
  AUX_COULOMB = 1,
  MINBAS = 2,
  HUECKEL = 3,
  IAO_LOCALIZATION = 4,
  SCF_DENS_GUESS = 5,
  AUX_CORREL = 6,
  ATOMIC_CHOLESKY = 7,
  ATOMIC_COMPACT_CHOLESKY = 8,
  ERF_ATOMIC_CHOLESKY = 9,
  ERF_ATOMIC_COMPACT_CHOLESKY = 10,
  AUX_JK = 11
};
template<>
void resolve<BASIS_PURPOSES>(std::string& value, BASIS_PURPOSES& field);

/**************************************************************************************************/
/*                                           Basis                                                */
/**************************************************************************************************/
/**
 * For what a certain basis set should be used.\n
 * COULOMB:          Auxiliary basis used to calculate Coulomb contributions.
 * EXCHANGE:         Auxiliary basis used to calculate Exchange (and Coulomb) contributions.
 * LREXCHANGE:       Auxiliary basis used to calculate Long-range Exchange contributions.
 * CORRELATION       Auxiliary basis used to calculate correlation contributions (e.g. MP2)
 */
enum class AUX_BASIS_PURPOSES { COULOMB = 0, EXCHANGE = 1, LREXCHANGE = 2, CORRELATION = 3 };
template<>
void resolve<AUX_BASIS_PURPOSES>(std::string& value, AUX_BASIS_PURPOSES& field);

/**************************************************************************************************/
/*                                     Density Fitting                                            */
/**************************************************************************************************/
/**
 * How the electron repulsion integrals should be calculated.
 * NONE: The full ERI is calculated and evaluated.
 * RI: The resolution of the identity approach is used to approximate the ERIs or corresponding contributions.
 * CD: A Cholesky decomposition of the complete ERIs is to evaluate them.
 * ACD: The atomic Cholesky decomposition approach is used to approximate the ERIs or corresponding contributions.
 * ACCD: The atomic-compact Cholesky decomposition approach is used to approximate the ERIs or corresponding
 * contributions.
 *
 * (Old explanation:
 * In case of a split Coulomb and exchange part in the Fock update (e.g. in DFT no Exchange is
 * needed) these are the possibilities to evaluate the Coulomb part. In principle all complete
 * HARTREE_FOCK_POTENTIAL_CALCULATORS are possible, but also the famous RI approximation, which is a
 * way to reduce the scaling behavior of DFT from O(n^4) to O(n^3).)
 */
enum class DENS_FITS { RI = 0, NONE = 1, CD = 2, ACD = 3, ACCD = 4 };
template<>
void resolve<DENS_FITS>(std::string& value, DENS_FITS& field);

/**************************************************************************************************/
/*                                     Extend Spherical ACD Shells                                            */
/**************************************************************************************************/
/**
 * In a pure spherical basis the use of a combined shell as the product of the base shells is an approximation
 * (If the combined shell is treated as a spherical shell). This can be circumvented by extending the
 * combined shell:
 *  NONE: no extension is applied
 *  SIMPLE: only one additional basis function is added for lower am.
 *  FIRST: Mixture of the COMPLETE (for the first correction) and SIMPLE procedure.
 *  COMPLETE: three additional basis functions are added for each basis functions of higher am.  (This
 *  is an experimental setting and can produce numeric problems for some systems!)
 */
enum class EXTEND_ACD { NONE = 0, SIMPLE = 1, FIRST = 2, COMPLETE = 3 };
template<>
void resolve<EXTEND_ACD>(std::string& value, EXTEND_ACD& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_BASISOPTIONS_H_ */
