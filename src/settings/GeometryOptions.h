/**
 * @file GeometryOptions.h
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

#ifndef SETTINGS_GEOMETRYOPTIONS_H_
#define SETTINGS_GEOMETRYOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                         Optimization                                           */
/**************************************************************************************************/
/**
 * Which algorithm to use for the geometry optimization.\n
 * SD: Steepest descent.\n
 * BFGS: Broyden-Fletcher-Goldfarb-Shanno.\n
 */
enum class OPTIMIZATION_ALGORITHMS { SQNM = 0, BFGS = 1 };
template<>
void resolve<OPTIMIZATION_ALGORITHMS>(std::string& value, OPTIMIZATION_ALGORITHMS& field);
/**************************************************************************************************/
/*                                    Geometry optimization                                       */
/**************************************************************************************************/
/**
 * How to calculate geometry gradients.\n
 * NUMERICAL: Calculate them numerically with finite steps.\n
 * ANALYTICAL: Calculate them analytically.
 */
enum class GRADIENT_TYPES { NUMERICAL = 0, ANALYTICAL = 1 };
template<>
void resolve<GRADIENT_TYPES>(std::string& value, GRADIENT_TYPES& field);

/**
 * Which algorithm to use for the geometry optimization.\n
 * GROUNDSTATE: Ground state optimization.\n
 * TS:          Transition state optimization
 */
enum class GEOMETRY_OPTIMIZATION_TYPES { GROUNDSTATE = 0, TS = 1 };
template<>
void resolve<GEOMETRY_OPTIMIZATION_TYPES>(std::string& value, GEOMETRY_OPTIMIZATION_TYPES& field);
/**
 * How to calculate the Hessian.\n
 * NUMERICAL: Calculate numerically with finite steps.\n
 * ANALYTICAL: Calculate analytically.
 */
enum class HESSIAN_TYPES { NUMERICAL = 0, ANALYTICAL = 1 };
template<>
void resolve<HESSIAN_TYPES>(std::string& value, HESSIAN_TYPES& field);

} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_GEOMETRYOPTIONS_H_ */
