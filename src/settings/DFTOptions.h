/**
 * @file DFTOptions.h
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef SETTINGS_DFTOPTIONS_H_
#define SETTINGS_DFTOPTIONS_H_

/* Include Serenity Internal Headers */
#include "dft/functionals/BasicFunctionals.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "settings/Options.h"

namespace Serenity {
namespace Options {

/**************************************************************************************************/
/*                                         Functionals                                            */
/**************************************************************************************************/
/* See CompositeFunctionals.h/.cpp for definitions and implementations */
template<>
void resolve<CompositeFunctionals::FUNCTIONALS>(std::string& value, CompositeFunctionals::FUNCTIONALS& field);
template<>
void resolve<CompositeFunctionals::XCFUNCTIONALS>(std::string& value, CompositeFunctionals::XCFUNCTIONALS& field);
template<>
void resolve<CompositeFunctionals::KINFUNCTIONALS>(std::string& value, CompositeFunctionals::KINFUNCTIONALS& field);
// Resolves the list inputs for XC, Kin functionals and embedding modes
template<>
void resolve<std::vector<CompositeFunctionals::XCFUNCTIONALS>>(std::string& value,
                                                               std::vector<CompositeFunctionals::XCFUNCTIONALS>& field);
template<>
void resolve<std::vector<CompositeFunctionals::KINFUNCTIONALS>>(std::string& value,
                                                                std::vector<CompositeFunctionals::KINFUNCTIONALS>& field);

template<>
void resolve<BasicFunctionals::BASIC_FUNCTIONALS>(std::string& value, BasicFunctionals::BASIC_FUNCTIONALS& field);

template<>
void resolve<CompositeFunctionals::IMPLEMENTATIONS>(std::string& value, CompositeFunctionals::IMPLEMENTATIONS& field);

template<>
void resolve<std::vector<BasicFunctionals::BASIC_FUNCTIONALS>>(std::string& value,
                                                               std::vector<BasicFunctionals::BASIC_FUNCTIONALS>& field);

/**************************************************************************************************/
/*                                         Dispersion                                             */
/**************************************************************************************************/
/**
 * Correction levels for the DFT-D corrections by Grimme et. al.
 * This includes all corrections (energy, gradient, hessian)
 *
 * NONE:     No correction.
 * D3:       The third set of parameters (with zero damping, also called D3(0)).
 * D3ABC:    The third set of parameters (with zero damping, also called D3(0)) and 3 center correction term.
 * D3BJ:     The third set of parameters with Becke-Johnson damping.
 * D3BJABC: The third set of parameters with Becke-Johnson damping and 3 center correction term.
 */
enum class DFT_DISPERSION_CORRECTIONS { NONE = 0, D3 = 1, D3ABC = 2, D3BJ = 3, D3BJABC = 4 };
template<>
void resolve<DFT_DISPERSION_CORRECTIONS>(std::string& value, DFT_DISPERSION_CORRECTIONS& field);

} /* namespace Options */
} /* namespace Serenity */
#endif /* SETTINGS_DFTOPTIONS_H_ */
