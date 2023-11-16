/**
 * @file GridOptions.h
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

#ifndef SETTINGS_GRIDOPTIONS_H_
#define SETTINGS_GRIDOPTIONS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"

namespace Serenity {
namespace Options {
/**************************************************************************************************/
/*                                            Grid                                                */
/**************************************************************************************************/
/**
 * How the radial part of the integration grid(s) is set up\n
 * EQUI: equidistant points
 * All others are schemes named after their authors,
 *  see the grid routines for the references.
 */
enum class RADIAL_GRID_TYPES {
  BECKE = 0,
  HANDY = 1,
  AHLRICHS = 2,
  KNOWLES = 3,
  EQUI = 4 //, MULTIEXP, EQUI
};
template<>
void resolve<RADIAL_GRID_TYPES>(std::string& value, RADIAL_GRID_TYPES& field);
/**
 * How the spherical part of the integration grid(s) is set up.
 *  All schemes are named after their authors.
 */
enum class SPHERICAL_GRID_TYPES { LEBEDEV = 0 };
template<>
void resolve<SPHERICAL_GRID_TYPES>(std::string& value, SPHERICAL_GRID_TYPES& field);
/**
 * How the atomic cells are determined - for each atom a grid is set up, they are then added
 * upon each other and recieve an additional weight prefactor based on the partitioning\n
 * BECKE:   the fuzzy cells (see A.D. Becke, J.Chem.Phys. 88 (1988), 2547.)\n
 * VORONOI: sharp cuts between the atoms (it is contrasted to the fuzzy cells in the paper above).\n
 * SSF:     grid according to Chem. Phys. Lett 257, (1996) 213-223
 */
enum class GRID_TYPES { BECKE = 0, VORONOI = 1, SSF = 2 };
template<>
void resolve<GRID_TYPES>(std::string& value, GRID_TYPES& field);
/**
 * For what a certain numerical grid should be used\n
 * DEFAULT: The standard integration grid to, e.g, calculate the DFT energy\n
 * SMALL:   A smaller integration grid for faster computation. This is typically used during an SCF
 *          in DFT calculations as long as the SCF is not converged\n
 * PLOT:    A grid to plot data on to for visualization purposes. Typically a cubical (equidistant)
 *          grid.
 */
enum class GRID_PURPOSES { DEFAULT = 0, SMALL = 1, PLOT = 2 };
template<>
void resolve<GRID_PURPOSES>(std::string& value, GRID_PURPOSES& field);
} /* namespace Options */
} /* namespace Serenity */

#endif /* SETTINGS_GRIDOPTIONS_H_ */
