/**
 * @file DFTOptions.cpp
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

#include "settings/DFTOptions.h"
#include "settings/Options.h"

namespace Serenity {
namespace Options {

template<>
void resolve<DFT_DISPERSION_CORRECTIONS>(std::string& value, DFT_DISPERSION_CORRECTIONS& field) {
  static const std::map<std::string, DFT_DISPERSION_CORRECTIONS> m = {{"NONE", DFT_DISPERSION_CORRECTIONS::NONE},
                                                                      {"D3", DFT_DISPERSION_CORRECTIONS::D3},
                                                                      {"D3ABC", DFT_DISPERSION_CORRECTIONS::D3ABC},
                                                                      {"D3BJ", DFT_DISPERSION_CORRECTIONS::D3BJ},
                                                                      {"D3BJABC", DFT_DISPERSION_CORRECTIONS::D3BJABC}};
  check(m, value, field);
}

} /* namespace Options */
} /* namespace Serenity */
