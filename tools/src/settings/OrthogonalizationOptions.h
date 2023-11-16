/**
 * @file OrthogonalizationOptions.h
 *
 * @author Anja Massolle
 * @date Aug 06, 2020
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

#ifndef SETTINGS_ORTHOGONALIZATIONOPTIONS_H_
#define SETTINGS_ORTHOGONALIZATIONOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                    Orthogonalization                                           */
/**************************************************************************************************/

/**
 * Orbital orthogonalization methods.
 *  NONE:                No orthogonalization criterion.
 */
enum class ORTHOGONALIZATION_ALGORITHMS { NONE = 0, LOEWDIN = 1, PIPEK = 2, BROER = 3 };
template<>
void resolve<ORTHOGONALIZATION_ALGORITHMS>(std::string& value, ORTHOGONALIZATION_ALGORITHMS& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_ORTHOGONALIZATIONOPTIONS_H_ */