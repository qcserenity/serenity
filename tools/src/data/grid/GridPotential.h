/**
 * @file   GridPotential.h
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   30. Juni 2015, last rework May 7, 2017 (JU)
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
#ifndef GRIDPOTENTIAL_H
#define GRIDPOTENTIAL_H
/* Include Serenity Internal Headers */
#include "data/grid/GridData.h"

namespace Serenity {
/* Forward Declarations */

/**
 * @class Serenity::GridPotential GridPotential.h
 * @brief Marker for scalar potentials which are represented on an integration grid. See GridData<SCFMode>
 */
template<Options::SCF_MODES SCFMode>
using GridPotential = GridData<SCFMode>;
} /* namespace Serenity */
#endif /* GRIDPOTENTIAL_H */
