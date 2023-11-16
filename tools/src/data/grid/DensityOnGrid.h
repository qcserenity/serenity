/**
 * @file   DensityOnGrid.h
 * @author Thomas Dresselhaus, Jan Unsleber
 *
 * @date   29. Dezember 2014, last rework May 7, 2017 (JU)
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
#ifndef DENSITYONGRID_H
#define DENSITYONGRID_H
/* Include Serenity Internal Headers */
#include "data/grid/GridData.h"

namespace Serenity {
/**
 * @class Serenity::DensityOnGrid DensityOnGrid.h
 * @brief Marker for the electron density represented on an integration grid. See GridData<SCFMode>
 */
template<Options::SCF_MODES SCFMode>
using DensityOnGrid = GridData<SCFMode>;
} /* namespace Serenity */
#endif /* DENSITYONGRID_H */
