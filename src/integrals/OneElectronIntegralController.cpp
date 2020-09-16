/**
 * @file   OneElectronIntegralController.cpp
 *
 * @date   30. Dezember 2013, 19:51
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
#include "integrals/OneElectronIntegralController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libecpint.h"
#include "integrals/wrappers/Libint.h"
//#include "data/matrices/MatrixInBasis.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

OneElectronIntegralController::OneElectronIntegralController(std::shared_ptr<BasisController> basisController,
                                                             std::shared_ptr<const Geometry> geometry)
  : _basisController(basisController), _geometry(geometry) {
  _basisController->addSensitiveObject(this->_self);
  // Check for ECPs
  _calcECPs = false;
  for (const auto& atom : _geometry->getAtoms()) {
    if (atom->usesECP())
      _calcECPs = true;
  }
}

OneElectronIntegralController::~OneElectronIntegralController() {
}

void OneElectronIntegralController::clearOneInts() {
  _overlapIntegrals.reset(nullptr);
  _oneElectronIntegrals.reset(nullptr);
  _ecpIntegrals.reset(nullptr);
  _oneElectronIntegralsTotal.reset(nullptr);
  _diplen.reset(nullptr);
  _dipvel.reset(nullptr);
  _angmom.reset(nullptr);
}

void OneElectronIntegralController::calcOverlapIntegrals() {
  takeTime("calculation of one-electron integrals");
  _overlapIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  auto& libint = Libint::getInstance();
  (*_overlapIntegrals) = libint.compute1eInts(libint2::Operator::overlap, _basisController);
  timeTaken(2, "calculation of one-electron integrals");
}

void OneElectronIntegralController::calcHCoreIntegrals() {
  takeTime("calculation of one-electron integrals");
  _oneElectronIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _oneElectronIntegralsTotal.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _ecpIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));

  auto& libint = Libint::getInstance();
  (*_oneElectronIntegrals) = libint.compute1eInts(libint2::Operator::kinetic, _basisController);
  (*_oneElectronIntegrals) += libint.compute1eInts(libint2::Operator::nuclear, _basisController, _geometry->getAtoms());
  _ecpIntegrals->setZero();
  if (_calcECPs)
    (*_ecpIntegrals) = Libecpint::computeECPIntegrals(_basisController, _geometry->getAtoms());
  *_oneElectronIntegralsTotal = *_ecpIntegrals + *_oneElectronIntegrals;
  timeTaken(2, "calculation of one-electron integrals");
}

void OneElectronIntegralController::calcDipoleLengths(Point gaugeOrigin) {
  takeTime("electric dipole integrals (length)");
  _diplen = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::emultipole1, 0, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(libint2::Operator::emultipole1, 0, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Matrix ints has 4 columns:
            // (i|1|j) -> (:,0),  (i|x|j) -> (:,1),   (i|y|j) -> (:,2),   (i|z|j) -> (:,3)
            // Electric dipole integrals in length rep (-r)
            (*_diplen)[0](mu, nu) = -ints((nj * k + l), 1);
            (*_diplen)[1](mu, nu) = -ints((nj * k + l), 2);
            (*_diplen)[2](mu, nu) = -ints((nj * k + l), 3);
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(libint2::Operator::emultipole1, 0, 2);
  timeTaken(2, "electric dipole integrals (length)");
}

void OneElectronIntegralController::calcDipoleVelocities(Point gaugeOrigin) {
  takeTime("electric dipole integrals (velocity)");
  _dipvel = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::emultipole1, 1, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(libint2::Operator::emultipole1, 1, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Matrix intsDeriv has 24 columns (Note: derivative-level 1!)
            // (i|1|dx j) -> (:, 0),   (i|x|dx j) -> (:, 1),   (i|y|dx j) -> (:, 2),   (i|z|dx j) -> (:, 3),
            // (i|1|dy j) -> (:, 4),   (i|x|dy j) -> (:, 5),   (i|y|dy j) -> (:, 6),   (i|z|dy j) -> (:, 7),
            // (i|1|dz j) -> (:, 8),   (i|x|dz j) -> (:, 9),   (i|y|dz j) -> (:,10),   (i|z|dz j) -> (:,11),

            // (dx i|1|j) -> (:,12),   (dx i|x|j) -> (:,13),   (dx i|y|j) -> (:,14),   (dx i|z|j) -> (:,15),
            // (dy i|1|j) -> (:,16),   (dy i|x|j) -> (:,17),   (dy i|y|j) -> (:,18),   (dy i|z|j) -> (:,19),
            // (dz i|1|j) -> (:,20),   (dz i|x|j) -> (:,21),   (dz i|y|j) -> (:,22),   (dz i|z|j) -> (:,23).

            // Electric dipole integrals in velocity rep -(i|nabla|j)
            (*_dipvel)[0](mu, nu) = -ints((nj * k + l), 0);
            (*_dipvel)[1](mu, nu) = -ints((nj * k + l), 4);
            (*_dipvel)[2](mu, nu) = -ints((nj * k + l), 8);
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(libint2::Operator::emultipole1, 1, 2);
  timeTaken(2, "electric dipole integrals (velocity)");
}

void OneElectronIntegralController::calcDipoleMagnetics(Point gaugeOrigin) {
  takeTime("magnetic-dipole integrals");
  _angmom = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::emultipole1, 1, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(libint2::Operator::emultipole1, 1, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Magnetic dipole integrals 0.5 (i|r x nabla|j)
            (*_angmom)[0](mu, nu) = 0.5 * (ints((nj * k + l), 10) - ints((nj * k + l), 7));
            (*_angmom)[1](mu, nu) = 0.5 * (ints((nj * k + l), 3) - ints((nj * k + l), 9));
            (*_angmom)[2](mu, nu) = 0.5 * (ints((nj * k + l), 5) - ints((nj * k + l), 2));
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(libint2::Operator::emultipole1, 1, 2);
  timeTaken(2, "magnetic-dipole integrals");
}

} /* namespace Serenity */
