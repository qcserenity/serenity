/**
 * @file   XyzFileToGeometryConverter.cpp
 *
 * @date   Mar 19, 2013
 * @author Thomas Dresselhaus
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
/* Include Class Header*/
#include "geometry/XyzFileToGeometryConverter.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "geometry/AtomTypeFactory.h"
#include "geometry/Geometry.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <sstream>
#include <vector>

namespace Serenity {

XyzFileToGeometryConverter::XyzFileToGeometryConverter(std::string filePath) : GeometryFileReader(filePath) {
}

std::unique_ptr<Geometry> XyzFileToGeometryConverter::readGeometryFromLoadedFile(unsigned int index) {
  std::vector<std::shared_ptr<Atom>> atoms;
  std::string line;
  unsigned int numberOfAtoms;
  unsigned int currentIndex = 0;

  while (true) {
    if (_file) {
      // The first entry in the first line should be an integer: the number of atoms.
      getline(_file, line);
    }
    else {
      // No more molecule there. Thus index is higher than avaliable.
      return nullptr;
    }
    std::istringstream iss(line);
    numberOfAtoms = 0;
    iss >> numberOfAtoms;
    // Further check
    if (numberOfAtoms == 0) {
      return nullptr;
    }
    // We don't care about the second line. It's the title line.
    if (_file) {
      getline(_file, line);
    }
    else {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    if (currentIndex++ == index) {
      // We arrived at the requested molecule.
      break;
    }
    else {
      // Skip over this molecule
      for (unsigned int i = 0; i < numberOfAtoms; ++i) {
        if (_file) {
          getline(_file, line);
        }
        else {
          throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
        }
      }
    }
  }
  // Loop over the atoms and read them in.
  for (unsigned int i = 0; i < numberOfAtoms; ++i) {
    if (_file) {
      getline(_file, line);
    }
    else {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    std::istringstream iss(line);
    std::string element;
    iss >> element;
    if (element == "") {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    double x, y, z;
    if (iss.good()) {
      iss >> x;
      x *= ANGSTROM_TO_BOHR;
    }
    else {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    if (iss.good()) {
      iss >> y;
      y *= ANGSTROM_TO_BOHR;
    }
    else {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    if (iss.good()) {
      iss >> z;
      z *= ANGSTROM_TO_BOHR;
    }
    else {
      throw SerenityError((std::string) "Error while reading xyz file " + _filePath);
    }
    atoms.push_back(std::make_shared<Atom>(AtomTypeFactory::getAtomType(element), x, y, z));
  }
  return std::unique_ptr<Geometry>(new Geometry(atoms));
}

} /* namespace Serenity */
