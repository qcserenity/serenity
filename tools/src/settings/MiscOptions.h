/**
 * @file MiscOptions.h
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

#ifndef SETTINGS_MISCOPTIONS_H_
#define SETTINGS_MISCOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                 System Partitioning                                            */
/**************************************************************************************************/
/**
 * Algorithms used for system partitioning based on orbital localization and predefined atom selections.
 *
 * ENFORCE_CHARGES:      Select orbitals for the subsystem until spin and charge of the subsystem are met.
 *                       The order of the selection is based on population analysis.
 * BEST_MATCH:           Select orbitals for the subsystem based on their predominant localization.
 * POPULATION_THRESHOLD: Select orbitals for one specific, prioritized subsystem based on a
 *                       orbital-population threshold.
 */
enum class SYSTEM_SPLITTING_ALGORITHM {
  ENFORCE_CHARGES,
  BEST_MATCH,
  POPULATION_THRESHOLD,
  SPADE,
  SPADE_ENFORCE_CHARGES
};
template<>
void resolve<SYSTEM_SPLITTING_ALGORITHM>(std::string& value, SYSTEM_SPLITTING_ALGORITHM& field);
/**************************************************************************************************/
/*                                Basis Set Truncation Algorithms                                 */
/**************************************************************************************************/
/**
 * Possible basis-set truncation algorithms.
 *  NONE:                     No truncation.
 *  NET_POPULATION:           Employ Mulliken net population criterion.
 *  PRIMITIVE_NET_POPULATION  Employ Mulliken net population criterion and keep only the given fraction.
 */
enum class BASIS_SET_TRUNCATION_ALGORITHMS { NONE = 0, NET_POPULATION = 1, PRIMITIVE_NET_POPULATION = 2 };
template<>
void resolve<BASIS_SET_TRUNCATION_ALGORITHMS>(std::string& value, BASIS_SET_TRUNCATION_ALGORITHMS& field);
/**************************************************************************************************/
/*                                        Gauge origin                                            */
/**************************************************************************************************/
/**
 * Different choices in the gauge origin.
 *  COM:    Center of mass.
 *  ORIGIN: Origin of the cartesian-coordinate system.
 */
enum class GAUGE_ORIGIN { COM, ORIGIN };
template<>
void resolve<GAUGE_ORIGIN>(std::string& value, GAUGE_ORIGIN& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_MISCOPTIONS_H_ */
