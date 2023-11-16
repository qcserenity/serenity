/**
 * @file   DensityMatrix.h
 * @author Jan Unsleber, Thomas Dresselhaus
 * @date   rework on April 09. 2017
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
#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"

namespace Serenity {
/**
 * @class Serenity::DensityMatrix DensityMatrix.h
 * @brief Marker for the reduced one-particle density matrix. See SpinPolarizedMatrixInBasis.
 *
 * Also called bond-order-charge-density matrix. If you don't understand what is going on here,
 * also take a look into the header files for MatrixInBasis and for SpinPolarizedData.
 */

template<Options::SCF_MODES SCFMode>
using DensityMatrix = MatrixInBasis<SCFMode>;
} /* namespace Serenity */
#endif /* DENSITYMATRIX_H */
