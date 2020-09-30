/**
 * @file   MatrixInBasis.cpp
 * @author Jan Unsleber
 * @date   Sep 24, 2020
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

#include "data/matrices/MatrixInBasis.h"
#include "io/HDF5.h"

namespace Serenity {

void MatrixInBasis<UNRESTRICTED>::toHDF5(std::string fBaseName, std::string matrixName) {
  std::string name = fBaseName + ".mat.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, matrixName.c_str(), this->alpha);
  HDF5::save(file, matrixName.c_str(), this->beta);
  HDF5::save_scalar_attribute(file, "basisSetName", _basisController->getBasisString());
}

void MatrixInBasis<RESTRICTED>::toHDF5(std::string fBaseName, std::string matrixName) {
  std::string name = fBaseName + ".mat.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, matrixName.c_str(), *this);
  HDF5::save_scalar_attribute(file, "basisSetName", _basisController->getBasisString());
}

} /* namespace Serenity */
