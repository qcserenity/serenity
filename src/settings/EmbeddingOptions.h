/**
 * @file EmbeddingOptions.h
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

#ifndef SETTINGS_EMBEDDINGOPTIONS_H_
#define SETTINGS_EMBEDDINGOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

/**
 * The possible treatments of Pauli-Repulsion for subsystem based methods.
 *  NONE:                   None.
 *  NADD_FUNC:              Non-additive kinetic energy functional.
 *  LEVELSHIFT:             Level-shift of occupied environment orbitals.
 *  HUZINAGA:               Solve the Huzinaga equation.
 *  HOFFMANN:               Solve Hoffmann's equation.
 *  RECONSTRUCTION:         Reconstruct the effective potential.
 *  FERMI_SHIFTED_HUZINAGA: Solve a shifted Huzinaga equation.
 *  LOEWDIN                 Evaluate contributions arising from an expanded overlap matrix.
 *  ALMO:                   Absolutely localized molecular orbitals.
 */
enum class KIN_EMBEDDING_MODES {
  NONE = 0,
  NADD_FUNC = 1,
  LEVELSHIFT = 2,
  HUZINAGA = 3,
  HOFFMANN = 4,
  RECONSTRUCTION = 5,
  FERMI_SHIFTED_HUZINAGA = 6,
  LOEWDIN = 7,
  ALMO = 8
};
template<>
void resolve<KIN_EMBEDDING_MODES>(std::string& value, KIN_EMBEDDING_MODES& field);

template<>
void resolve<std::vector<KIN_EMBEDDING_MODES>>(std::string& value, std::vector<KIN_EMBEDDING_MODES>& field);
/**
 * Embedding modes which can be used for the BS-DFT task:
 * NONE: KS-DFT calculation
 * ISOLATED :
 * FDE : Frozen Density Embedding
 * FAT : Freeze and Thaw embedding
 */
enum class EMBEDDING_SCHEME { NONE = 0, ISOLATED = 1, FDE = 2, FAT = 3 };
template<>
void resolve<EMBEDDING_SCHEME>(std::string& value, EMBEDDING_SCHEME& field);

} /* namespace Options */
} /* namespace Serenity */
#endif /* SETTINGS_EMBEDDINGOPTIONS_H_ */
