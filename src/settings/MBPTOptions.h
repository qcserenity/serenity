/**
 * @file MBPTOptions.h
 *
 * @author Johannes Toelle
 * @date May 27, 2020
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef SETTINGS_MBPTOPTIONS_H_
#define SETTINGS_MBPTOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                           MBPT                                                */
/**************************************************************************************************/
/**
 * Possible types for Many-Body Perturbation Theory (MBPT) calculation to determine quasiparticle states:
 * GW : GW calculation
 * RPA: Random-Phase-Approximation
 */
enum class MBPT { GW = 0, RPA = 1 };
template<>
void resolve<MBPT>(std::string& value, MBPT& field);
/**************************************************************************************************/
/*                                           GW                                                */
/**************************************************************************************************/
/**
 * Possible types for GW calculation to determine quasiparticle states:
 * CD : Contour deformation
 * AC : Analytic Continuation
 * ANALYTIC: Full analytic evaluation
 */
enum class GWALGORITHM { CD = 0, AC = 1, ANALYTIC = 2 };
template<>
void resolve<GWALGORITHM>(std::string& value, GWALGORITHM& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_MBPTOPTIONS_H_ */