/**
 * @file LRSCFOptions.h
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

#ifndef SETTINGS_LRSCFOPTIONS_H_
#define SETTINGS_LRSCFOPTIONS_H_

#include "settings/Options.h"

namespace Serenity {
namespace Options {

/**************************************************************************************************/
/*                                           LRSCF                                                */
/**************************************************************************************************/
/**
 * Possible types for LRSCF problems to determine eigenvalue solving procedure:
 * TDA/CIS : AX = Xw (is Hermitian, uses Davidson)
 * TDDFT : sqrt(A-B)(A+B)sqrt(A-B) sqrt-(X+Y)
 *           = sqrt-(X+Y) w^2 (is Hermitian, uses Davidson)
 * RPA : (A+B)(X+Y) = (X-Y)w
 *       (A-B)(X-Y) = (X+Y)w (is non-Hermitian, uses a modifed OJJ)
 * The latter one includes TDHF, hybrid TDDFT, FDEc with
 * external orthogonality and supermolecular TDDFT with local orbitals.
 */
enum class RESPONSE_PROBLEM { TDA = 0, TDDFT = 1, RPA = 2 };
template<>
void resolve<RESPONSE_PROBLEM>(std::string& value, RESPONSE_PROBLEM& field);
/**
 * Type of the LRSCF calculation.
 * ISOLATED:  Only the fragment.
 * UNCOUPLED: Ground state embedding is considered.
 * COUPLED:   Full coupling of the fragments.
 */
enum class LRSCF_TYPE { ISOLATED = 0, UNCOUPLED = 1, COUPLED = 2 };
template<>
void resolve<LRSCF_TYPE>(std::string& value, LRSCF_TYPE& field);
/**
 * Options for integral evaluation.
 * NUMERICAL:  Numerical evaluation.
 * ANALYTICAL: Analytical evaluation.
 */
enum class INTEGRAL_TYPE { NUMERICAL = 0, ANALYTICAL = 1 };
template<>
void resolve<INTEGRAL_TYPE>(std::string& value, INTEGRAL_TYPE& field);
/**
 * Gauge options for dipole integrals.
 * LENGTH:   Length gauge.
 * VELOCITY: Velocity gauge.
 */
enum class GAUGE { LENGTH = 0, VELOCITY = 1 };
template<>
void resolve<GAUGE>(std::string& value, GAUGE& field);
/**
 * Multiplicity of excited state in LRSCF calculation
 * Singlet: 2*S + 1 = 1
 * Triplet: 2*S + 1 = 3
 */
enum class MULTIPLICITY { SINGLET = 0, TRIPLET = 1 };
template<>
void resolve<MULTIPLICITY>(std::string& value, MULTIPLICITY& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_LRSCFOPTIONS_H_ */
