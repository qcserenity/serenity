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

Atom::Atom(std::shared_ptr<const AtomType> atomType, const double x, const double y, const double z,
           const std::pair<std::string, std::vector<std::shared_ptr<Shell>>> basisFunctions)
  : Atom(atomType, x, y, z) {
  _primaryBasisLabel = basisFunctions.first;
  _associatedBasis.insert(basisFunctions);
}

Atom::Atom(std::shared_ptr<const AtomType> atomType, const double x, const double y, const double z)
  : Point(x, y, z), _atomType(atomType), _gradientsUpToDate(false), _gradient{0.0, 0.0, 0.0}, _primaryBasisLabel("-") {
}

Atom::Atom(const std::string symbol, const double x, const double y, const double z)
  : Atom(AtomTypeFactory::getAtomType(symbol), x, y, z) {
}

unsigned int Atom::getNBasisFunctions() const {
  return this->getNBasisFunctions(_primaryBasisLabel);
}

unsigned int Atom::getNBasisFunctions(const std::string label) const {
  if (_associatedBasis.size() > 0) {
    return _associatedBasis.at(label).size();
  }
  else {
    throw SerenityError((std::string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
}

bool Atom::operator==(Atom rhs) {
  bool same = true;
  auto lhsName = this->getAtomType()->getName();
  auto rhsName = rhs.getAtomType()->getName();
  same = same && this->_atomType->getElementSymbol() == rhs.getAtomType()->getElementSymbol();
  // Some structure editors save the coordinates only in single precision. Therefore, it makes
  // little sense to check for differences smaller than 5e-5. To be honest, even this threshold
  // may be too strict.
  same = same && isEqual(this->getX(), rhs.getX(), 5e-5);
  same = same && isEqual(this->getY(), rhs.getY(), 5e-5);
  same = same && isEqual(this->getZ(), rhs.getZ(), 5e-5);
  return same;
}

std::vector<std::shared_ptr<Shell>>& Atom::getBasisFunctions() {
  return this->getBasisFunctions(_primaryBasisLabel);
}

const std::vector<std::shared_ptr<Shell>>& Atom::getBasisFunctions() const {
  return this->getBasisFunctions(_primaryBasisLabel);
}

std::vector<std::shared_ptr<Shell>>& Atom::getBasisFunctions(const std::string label) {
  if (_associatedBasis.find(label) == _associatedBasis.end()) {
    throw SerenityError((std::string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
  else {
    return _associatedBasis.at(label);
  }
}

bool Atom::basisFunctionsExist(std::string label) {
  return !(_associatedBasis.find(label) == _associatedBasis.end());
}

const std::vector<std::shared_ptr<Shell>>& Atom::getBasisFunctions(const std::string label) const {
  if (_associatedBasis.size() > 0) {
    return _associatedBasis.at(label);
  }
  else {
    throw SerenityError((std::string) "You are requesting a non-existent basis from an atom. Label: " + label);
  }
}

void Atom::addGrid(const std::pair<std::string, AtomGrid*> newGrid, const bool isPrimary) {
  _grids.insert(newGrid);
  if (isPrimary) {
    _primaryGridLabel = newGrid.first;
  }
}

AtomGrid* Atom::getGrid() {
  return getGrid(_primaryGridLabel);
}

AtomGrid* Atom::getGrid(const std::string label) {
  return _grids.at(label);
}

void Atom::setX(const double x) {
  this->at(0) = x;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setX(x);
  _gradientsUpToDate = false;
  if (_corePotential) {
    _corePotential->setPos(this->getX(), this->getY(), this->getZ());
  }
  notifyObjects();
}

void Atom::setY(const double y) {
  this->at(1) = y;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setY(y);
  _gradientsUpToDate = false;
  if (_corePotential) {
    _corePotential->setPos(this->getX(), this->getY(), this->getZ());
  }
  notifyObjects();
}

void Atom::setZ(const double z) {
  this->at(2) = z;
  for (auto& basis : _associatedBasis)
    for (auto& basisFunction : basis.second)
      basisFunction->setZ(z);
  _gradientsUpToDate = false;
  if (_corePotential) {
    _corePotential->setPos(this->getX(), this->getY(), this->getZ());
  }
  notifyObjects();
}

void Atom::addToX(double add_to_x) {
  setX(this->at(0) + add_to_x);
}

void Atom::addToY(double add_to_y) {
  setY(this->at(1) + add_to_y);
}

void Atom::addToZ(double add_to_z) {
  setZ(this->at(2) + add_to_z);
}

std::shared_ptr<libecpint::ECP> Atom::getCorePotential() const {
  assert(_primaryBasisLabel != "-");
  assert(_corePotential);
  return _corePotential;
}
void Atom::addBasis(std::pair<std::string, std::vector<std::shared_ptr<Shell>>> newBasis,
                    std::shared_ptr<libecpint::ECP> ecp, unsigned int nECPElectrons) {
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
  if (_nECPElectrons == 0)
    _nECPElectrons = nECPElectrons;
  // TODO it may look logical to execute the notifyObjects() here, however,
  // this is currently not possible, since this would happen during basis
  // construction. At the current moment Serenity expects that the basis
  // set is built in one go without any notifys deleting it mid construction.
  // notifyObjects();
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
