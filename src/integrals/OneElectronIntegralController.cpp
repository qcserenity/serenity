/**
 * @file   OneElectronIntegralController.cpp
 * 
 * @date   30. Dezember 2013, 19:51
 * @author Thomas Dresselhaus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "integrals/OneElectronIntegralController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "integrals/wrappers/Libecpint.h"
//#include "data/matrices/MatrixInBasis.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>


namespace Serenity {
using namespace std;

OneElectronIntegralController::OneElectronIntegralController(
    shared_ptr<BasisController> basisController,
    shared_ptr<const Geometry> geometry) :
      _basisController(basisController),
      _geometry(geometry) {
  _basisController->addSensitiveObject(this->_self);
  // Check for ECPs
  _calcECPs = false;
  for (const auto& atom : _geometry->getAtoms()) {
    if (atom->usesECP()) _calcECPs = true;
  }
}

OneElectronIntegralController::~OneElectronIntegralController(){
}

void OneElectronIntegralController::clearOneInts() {
  _overlapIntegrals.reset(nullptr);
  _oneElectronIntegrals.reset(nullptr);
  _ecpIntegrals.reset(nullptr);
  _oneElectronIntegralsTotal.reset(nullptr);
}

double OneElectronIntegralController::S(unsigned int i, unsigned int j) {
  if (!_overlapIntegrals) {
    calcOverlapIntegrals();
  }
  return (*_overlapIntegrals)(i, j);
}

double OneElectronIntegralController::h(unsigned int i, unsigned int j) {
  if (!_oneElectronIntegralsTotal) {
    calcHCoreIntegrals();
  }
  return (*_oneElectronIntegralsTotal)(i,j);
}


void OneElectronIntegralController::calcOverlapIntegrals() {
  takeTime("calculation of one-electron integrals");
  _overlapIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  auto& libint = Libint::getInstance();
  (*_overlapIntegrals) = libint.compute1eInts(libint2::Operator::overlap, _basisController);
  timeTaken(2,"calculation of one-electron integrals");
}

void OneElectronIntegralController::calcHCoreIntegrals() {
  takeTime("calculation of one-electron integrals");
  _oneElectronIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _oneElectronIntegralsTotal.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _ecpIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));

  auto& libint = Libint::getInstance();
  (*_oneElectronIntegrals) = libint.compute1eInts(libint2::Operator::kinetic, _basisController);
  (*_oneElectronIntegrals) += libint.compute1eInts(libint2::Operator::nuclear, _basisController,_geometry->getAtoms());
  _ecpIntegrals->setZero();
  if (_calcECPs) (*_ecpIntegrals) = Libecpint::computeECPIntegrals(_basisController,_geometry->getAtoms());
  *_oneElectronIntegralsTotal = *_ecpIntegrals + *_oneElectronIntegrals;
  timeTaken(2,"calculation of one-electron integrals");
}

void OneElectronIntegralController::printOneElectronIntegrals() {
  if (!_oneElectronIntegralsTotal) {
    std::cout << "Cannot print one-electron integrals.\n"
         "None have yet been calculated for this class."
         << std::endl;
  } else {
    std::cout << *_oneElectronIntegralsTotal << std::endl;
  }
}
void OneElectronIntegralController::printOverlapIntegrals() {
  if (!_overlapIntegrals) {
    std::cout << "Cannot print overlap integrals.\n"
         "None have yet been calculated for this class."
         << std::endl;
  } else {
    std::cout << *_overlapIntegrals << std::endl;
  }
}

} /* namespace Serenity */
