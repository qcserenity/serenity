/**
 * @file ElectronicStructureOptions.h
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

#ifndef SETTINGS_ELECTRONICSTRUCTUREOPTIONS_H_
#define SETTINGS_ELECTRONICSTRUCTUREOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**
 * The file formats supported for the OrbitalsIO task.
 *   TURBOMOLE: Turbomole ASCII-MOS.
 *   SERENITY:  Serenity HDF5 files.
 *   MOLPRO:    Molpro-xml orbital file.
 *   MOLCAS:    Molcas hdf5 file format (.scf.h5 files).
 *   MOLDEN:    Molden file format.
 */
enum class ORBITAL_FILE_TYPES { SERENITY = 0, TURBOMOLE = 1, MOLPRO = 2, MOLCAS = 3, MOLDEN = 4 };
template<>
void resolve<ORBITAL_FILE_TYPES>(std::string& value, ORBITAL_FILE_TYPES& field);
/**
 * The file formats supported for the WriteIntegrals task.
 *   ASCII: ASCII files.
 *   HDF5:  HDF5 files.
 */
enum class INTEGRAL_FILE_TYPES { ASCII = 0, HDF5 = 1 };
template<>
void resolve<INTEGRAL_FILE_TYPES>(std::string& value, INTEGRAL_FILE_TYPES& field);
/**
 * The type of SCF calculation that is made\n
 * RESTRICTED: all electrons are paired, works only for even numbers of electrons\n
 * UNRESTRICTED: uneven numbers of electrons are allowed
 */
enum class SCF_MODES { RESTRICTED = 0, UNRESTRICTED = 1 };
template<>
void resolve<SCF_MODES>(std::string& value, SCF_MODES& field);
/**
 * The type of ROHF method (both constrained UHF variants.
 * NONE: No ROHF.
 * CUHF: J. Chem. Phys. 133, 141102 (2010).
 * SUHF: Chem. Phys. Lett. 183, 423 (1991).
 */
enum class ROHF_TYPES { NONE = 0, CUHF = 1, SUHF = 2 };
template<>
void resolve<ROHF_TYPES>(std::string& value, ROHF_TYPES& field);
/**
 * The electronic structure theory which is used on the single particle level.
 * HF:  Hartree-Fock
 * DFT: Density-functional theory.
 */
enum class ELECTRONIC_STRUCTURE_THEORIES { HF = 0, DFT = 1 };
template<>
void resolve<ELECTRONIC_STRUCTURE_THEORIES>(std::string& value, ELECTRONIC_STRUCTURE_THEORIES& field);
/**
 * Print level settings. The print level names should make any further description unnecessary.
 *  MINIMUM
 *  NORMAL
 *  VERBOSE
 *  DEBUG
 */
enum class GLOBAL_PRINT_LEVELS { MINIMUM, NORMAL, VERBOSE, DEBUGGING };
template<>
void resolve<GLOBAL_PRINT_LEVELS>(std::string& value, GLOBAL_PRINT_LEVELS& field);

} /* namespace Options */
static constexpr Options::SCF_MODES RESTRICTED = Options::SCF_MODES::RESTRICTED;
static constexpr Options::SCF_MODES UNRESTRICTED = Options::SCF_MODES::UNRESTRICTED;
} /* namespace Serenity */
#endif /* SETTINGS_ELECTRONICSTRUCTUREOPTIONS_H_ */
