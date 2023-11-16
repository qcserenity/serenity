/**
 * @file   FockMatrix.h
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
#ifndef FOCKMATRIX_H
#define FOCKMATRIX_H
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"

namespace Serenity {
/**
 * @class Serenity::FockMatrix FockMatrix.h
 * @brief Marker for the matrix representation of the operator used in one-particle equations.
 *
 * Usage: see SpinPolarizedMatrixInBasis.\n
 * The many-body problem (i.e. many electrons) we want to solve in electronic structure theory is
 * too complicated to be solved directly. Instead, at least in HartreeFock and Kohn--Sham Dft, an
 * effective potential is formed based on a guessed electronic structure and used in one-particle
 * equations, the Hartree--Fock equations (or Kohn--Sham equations). The operator used in these
 * equations is called the Fock operator and containins all information about the effective
 * potential, other potential terms and also the kinetic energy operator. This class is the matrix
 * form of that operator.\n
 *
 * Also take a look into the header file for MatrixInBasis.
 */
template<Options::SCF_MODES SCFMode>
using FockMatrix = MatrixInBasis<SCFMode>;
} /* namespace Serenity */
#endif /* FOCKMATRIX_H */
