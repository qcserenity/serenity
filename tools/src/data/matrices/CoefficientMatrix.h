/**
 * @file   CoefficientMatrix.h
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
#ifndef COEFFICIENTMATRIX_H
#define COEFFICIENTMATRIX_H
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"

namespace Serenity {
/**
 * @class Serenity::CoefficientMatrix CoefficientMatrix.h
 * @brief Marker for matrices containing data about the composition of orbitals
 *
 * A set of column vectors determining the orbitals in connection with a basis. The values inside
 * the vectors specify how much a basis function contributes to the orbital.
 * Use as coefficientMatrix(<basis function>,<orbital>).
 *
 * For a more efficient use you can also access the underlying data directly as is usual for
 * matrices.
 */
template<Options::SCF_MODES SCFMode>
using CoefficientMatrix = MatrixInBasis<SCFMode>;
} /* namespace Serenity */
#endif /* COEFFICIENTMATRIX_H */
