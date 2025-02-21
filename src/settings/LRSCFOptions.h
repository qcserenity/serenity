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

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {

/**************************************************************************************************/
/*                                           LRSCF                                                */
/**************************************************************************************************/
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
 * Excited-state wavefunction model.
 */
enum class LR_METHOD { TDA = 0, TDDFT = 1, CC2 = 2, CISDINF = 3, CISD = 4, ADC2 = 5 };
template<>
void resolve<LR_METHOD>(std::string& value, LR_METHOD& field);

/**
 * Different stability analyses:
 *
 * Real RHF -> Real RHF    : (A+B), singlet
 * Real RHF -> Real UHF    : (A+B), triplet
 * Real RHF -> Complex RHF : (A-B), singlet
 *
 * Real UHF -> Real UHF    : (A+B)
 * Real UHF -> Complex UHF : (A-B)
 *
 * Spin-Flip:
 *   1. Start from triplet reference: nElectrons(alpha) = nElectrons(beta) + 2
 *   2. Occupied reference: alpha
 *   3. Virtual reference: beta
 *   4. No coupling matrix contributions from Coulomb and XC Kernel (only HF exchange).
 *   5. Automatically done within the TDA.
 */
enum class STABILITY_ANALYSIS { NONE = 0, REAL = 1, NONREAL = 2, SPINFLIP = 3 };
template<>
void resolve<STABILITY_ANALYSIS>(std::string& value, STABILITY_ANALYSIS& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_LRSCFOPTIONS_H_ */
