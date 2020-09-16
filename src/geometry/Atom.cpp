/**
 * @file   Atom.cpp
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
#include "geometry/Atom.h"
/* Include Serenity Internal Headers */
#include "basis/Shell.h"
#include "geometry/AtomTypeFactory.h"
#include "grid/construction/AtomGrid.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <algorithm>
#include <libecpint/ecp.hpp>
#include <stdexcept>

namespace Serenity {
using namespace std;

Atom::Atom(shared_ptr<const AtomType> atomType, const double x, const double y, const double z,
           const pair<string, vector<shared_ptr<Shell>>> basisFunctions)
  : Atom(atomType, x, y, z) {
  _primaryBasisLabel = basisFunctions.first;
  _associatedBasis.insert(basisFunctions);
}

Atom::Atom(shared_ptr<const AtomType> atomType, const double x, const double y, const double z)
  : Point(x, y, z), _atomType(atomType), _gradientsUpToDate(false), _gradient{0.0, 0.0, 0.0}, _primaryBasisLabel("-") {
}

Atom::Atom(const std::string symbol, const double x, const double y, const double z)
  : Atom(AtomTypeFactory::getAtomType(symbol), x, y, z) {
}

unsigned int Atom::getNBasisFunctions() const {
  return this->getNBasisFunctions(_primaryBasisLabel);
}

unsigned int Atom::getNBasisFunctions(const string label) const {
  if (_associatedBasis.size() > 0) {
    return _associatedBasis.at(label).size();
  }
  else {
    throw SerenityError((string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
}

vector<shared_ptr<Shell>>& Atom::getBasisFunctions() {
  return this->getBasisFunctions(_primaryBasisLabel);
}

const vector<shared_ptr<Shell>>& Atom::getBasisFunctions() const {
  return this->getBasisFunctions(_primaryBasisLabel);
}

vector<shared_ptr<Shell>>& Atom::getBasisFunctions(const string label) {
  if (_associatedBasis.find(label) == _associatedBasis.end()) {
    throw SerenityError((string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
  else {
    return _associatedBasis.at(label);
  }
}

bool Atom::basisFunctionsExist(std::string label) {
  return !(_associatedBasis.find(label) == _associatedBasis.end());
}

const vector<shared_ptr<Shell>>& Atom::getBasisFunctions(const string label) const {
  if (_associatedBasis.size() > 0) {
    return _associatedBasis.at(label);
  }
  else {
    throw SerenityError((string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
}

void Atom::addGrid(const pair<string, AtomGrid*> newGrid, const bool isPrimary) {
  _grids.insert(newGrid);
  if (isPrimary) {
    _primaryGridLabel = newGrid.first;
  }
}

AtomGrid* Atom::getGrid() {
  return getGrid(_primaryGridLabel);
}

AtomGrid* Atom::getGrid(const string label) {
  return _grids.at(label);
}

void Atom::setX(const double x) {
  _x = x;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setX(x);
  _gradientsUpToDate = false;
  notifyObjects();
}

void Atom::setY(const double y) {
  _y = y;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setY(y);
  _gradientsUpToDate = false;
  notifyObjects();
}

void Atom::setZ(const double z) {
  _z = z;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setZ(z);
  _gradientsUpToDate = false;
  notifyObjects();
}

void Atom::addToX(double add_to_x) {
  setX(_x + add_to_x);
}

void Atom::addToY(double add_to_y) {
  setY(_y + add_to_y);
}

void Atom::addToZ(double add_to_z) {
  setZ(_z + add_to_z);
}

std::shared_ptr<libecpint::ECP> Atom::getCorePotential() const {
  assert(_primaryBasisLabel != "-");
  assert(_corePotential);
  return _corePotential;
}
void Atom::addBasis(std::pair<std::string, std::vector<std::shared_ptr<Shell>>> newBasis,
                    std::shared_ptr<libecpint::ECP> ecp, unsigned int nCoreElectrons) {
  // Make sure no primary basis has been specified before.
  if (_primaryBasisLabel != "-" and ecp->N != 0)
    throw SerenityError("Error: A set of effective core potential functions is supposed to be added to an atom which"
                        " already has a 'primary' basis. This is not allowed, because using ECPs changes the"
                        " effective charge and number of electrons. This may lead to errors if data is used which"
                        " has been calculated before without using ECPs.");
  this->addBasis(newBasis, true);
  if (ecp->N != 0) {
    assert(!_corePotential);
    _corePotential = ecp;
  }
  // Only change the number of core electrons if there was no ECP before.
  if (_nCoreElectrons == 0)
    _nCoreElectrons = nCoreElectrons;
  notifyObjects();
}

const Eigen::Vector3d& Atom::getGradient() const {
  if (_gradientsUpToDate) {
    return _gradient;
  }
  else {
    throw SerenityError("Gradients requested from an atom, but they were out of date.");
  }
}

} /* namespace Serenity */
