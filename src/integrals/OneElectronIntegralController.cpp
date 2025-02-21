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
#include "basis/Basis.h" //Loop shells.
#include "basis/BasisController.h"
#include "data/ExternalChargeController.h"
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libecpint.h"
#include "integrals/wrappers/Libint.h"
#include "math/linearAlgebra/MatrixFunctions.h" //symmetrize()
#include "misc/Timing.h"                        //Timings
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

OneElectronIntegralController::OneElectronIntegralController(std::shared_ptr<BasisController> basisController,
                                                             std::shared_ptr<const Geometry> geometry,
                                                             std::shared_ptr<ExternalChargeController> externalChargeController)
  : _basisController(basisController), _geometry(geometry), _externalChargeController(externalChargeController) {
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
  _kinIntegrals.reset(nullptr);
  _nucIntegrals.reset(nullptr);
  _diplen.reset(nullptr);
  _dipvel.reset(nullptr);
  _angmom.reset(nullptr);
  _extChargeIntegrals.reset(nullptr);
}

void OneElectronIntegralController::calcOverlapIntegrals() {
  takeTime("calculation of one-electron integrals");
  _overlapIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  auto& libint = Libint::getInstance();
  //  (*_overlapIntegrals) = symmetrize(libint.compute1eInts(LIBINT_OPERATOR::overlap, _basisController));
  (*_overlapIntegrals) = libint.compute1eInts(LIBINT_OPERATOR::overlap, _basisController);
  timeTaken(2, "calculation of one-electron integrals");
}

void OneElectronIntegralController::calcKinIntegrals() {
  takeTime("calculation of one-electron integrals");
  _kinIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  auto& libint = Libint::getInstance();
  *_kinIntegrals = libint.compute1eInts(LIBINT_OPERATOR::kinetic, _basisController);
  timeTaken(2, "calculation of one-electron integrals");
}
void OneElectronIntegralController::calcNucIntegrals() {
  takeTime("calculation of one-electron integrals");
  _nucIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  auto& libint = Libint::getInstance();
  *_nucIntegrals = libint.compute1eInts(LIBINT_OPERATOR::nuclear, _basisController, _geometry->getAtoms());
  timeTaken(2, "calculation of one-electron integrals");
}

void OneElectronIntegralController::calcECPIntegrals() {
  takeTime("calculation of one-electron integrals");
  _ecpIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _ecpIntegrals->setZero();
  if (_calcECPs)
    (*_ecpIntegrals) = symmetrize(Libecpint::computeECPIntegrals(_basisController, _geometry->getAtoms()));
  timeTaken(2, "calculation of one-electron integrals");
}

void OneElectronIntegralController::calcHCoreIntegrals() {
  _oneElectronIntegrals.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  _oneElectronIntegralsTotal.reset(new MatrixInBasis<RESTRICTED>(_basisController));
  (*_oneElectronIntegrals) = this->getKinIntegrals();
  (*_oneElectronIntegrals) += this->getNucIntegrals();
  *_oneElectronIntegralsTotal = this->getECPIntegrals() + *_oneElectronIntegrals;
  *_oneElectronIntegralsTotal = symmetrize(*_oneElectronIntegralsTotal);
  (*_oneElectronIntegrals) = symmetrize(*_oneElectronIntegrals);
}

void OneElectronIntegralController::calcDipoleLengths(Point gaugeOrigin) {
  takeTime("electric dipole integrals (length)");
  _diplen = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::emultipole1, 0, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(LIBINT_OPERATOR::emultipole1, 0, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Matrix ints has 4 columns:
            // (i|1|j) -> (:,0),  (i|x|j) -> (:,1),   (i|y|j) -> (:,2),   (i|z|j) -> (:,3)
            // Electric dipole integrals in length rep
            (*_diplen)[0](mu, nu) = -ints((nj * k + l), 1);
            (*_diplen)[1](mu, nu) = -ints((nj * k + l), 2);
            (*_diplen)[2](mu, nu) = -ints((nj * k + l), 3);
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(LIBINT_OPERATOR::emultipole1, 0, 2);
  timeTaken(2, "electric dipole integrals (length)");
}

void OneElectronIntegralController::calcDipoleVelocities(Point gaugeOrigin) {
  takeTime("electric dipole integrals (velocity)");
  _dipvel = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::emultipole1, 1, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(LIBINT_OPERATOR::emultipole1, 1, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Matrix intsDeriv has 24 columns (Note: derivative-level 1!):

            // (dx i|1|j) -> (:, 0),
            // (dx i|x|j) -> (:, 1),
            // (dx i|y|j) -> (:, 2),
            // (dx i|z|j) -> (:, 3),

            // (dy i|1|j) -> (:, 4),
            // (dy i|x|j) -> (:, 5),
            // (dy i|y|j) -> (:, 6),
            // (dy i|z|j) -> (:, 7),

            // (dz i|1|j) -> (:, 8),
            // (dz i|x|j) -> (:, 9),
            // (dz i|y|j) -> (:,10),
            // (dz i|z|j) -> (:,11),

            // (i|1|dx j) -> (:,12),
            // (i|x|dx j) -> (:,13),
            // (i|y|dx j) -> (:,14),
            // (i|z|dx j) -> (:,15),

            // (i|1|dy j) -> (:,16),
            // (i|x|dy j) -> (:,17),
            // (i|y|dy j) -> (:,18),
            // (i|z|dy j) -> (:,19),

            // (i|1|dz j) -> (:,20),
            // (i|x|dz j) -> (:,21),
            // (i|y|dz j) -> (:,22),
            // (i|z|dz j) -> (:,23).

            // Electric dipole integrals in velocity rep
            (*_dipvel)[0](mu, nu) = -ints((nj * k + l), 12);
            (*_dipvel)[1](mu, nu) = -ints((nj * k + l), 16);
            (*_dipvel)[2](mu, nu) = -ints((nj * k + l), 20);
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(LIBINT_OPERATOR::emultipole1, 1, 2);
  timeTaken(2, "electric dipole integrals (velocity)");
}

void OneElectronIntegralController::calcDipoleMagnetics(Point gaugeOrigin) {
  takeTime("magnetic-dipole integrals");
  _angmom = std::make_unique<std::vector<MatrixInBasis<RESTRICTED>>>(3, _basisController);

  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::emultipole1, 1, 2, gaugeOrigin);
  auto basis = _basisController->getBasis();
  for (unsigned i = 0; i < basis.size(); i++) {
    for (unsigned j = 0; j < basis.size(); j++) {
      Eigen::MatrixXd ints;
      if (libint.compute(LIBINT_OPERATOR::emultipole1, 1, *basis[i], *basis[j], ints)) {
        for (unsigned k = 0; k < basis[i]->getNContracted(); k++) {
          auto mu = _basisController->extendedIndex(i) + k;
          for (unsigned l = 0; l < basis[j]->getNContracted(); l++) {
            auto nu = _basisController->extendedIndex(j) + l;
            unsigned nj = basis[j]->getNContracted();
            // Magnetic dipole integrals
            (*_angmom)[0](mu, nu) = -0.5 * (ints((nj * k + l), 22) - ints((nj * k + l), 19));
            (*_angmom)[1](mu, nu) = -0.5 * (ints((nj * k + l), 15) - ints((nj * k + l), 21));
            (*_angmom)[2](mu, nu) = -0.5 * (ints((nj * k + l), 17) - ints((nj * k + l), 14));
          }
        }
      }
    } /* shell j */
  }   /* shell i */
  libint.finalize(LIBINT_OPERATOR::emultipole1, 1, 2);
  timeTaken(2, "magnetic-dipole integrals");
}
void OneElectronIntegralController::calcExternalChargeIntegrals() {
  if (!this->getExternalCharges().empty()) {
    auto& libint = Libint::getInstance();
    takeTime("DOI - analytical");
    auto shellPairData = _basisController->getDOIPrescreeningFactors();
    timeTaken(2, "DOI - analytical");
    _extChargeIntegrals = std::make_unique<MatrixInBasis<RESTRICTED>>(this->_basisController);
    *_extChargeIntegrals =
        libint.compute1eInts(LIBINT_OPERATOR::nuclear, this->_basisController, this->getExternalCharges(),
                             std::numeric_limits<double>::epsilon(), 10, Libint::getNPrimMax(), shellPairData);
    return;
  }
  _extChargeIntegrals = std::make_unique<MatrixInBasis<RESTRICTED>>(this->_basisController);
  _extChargeIntegrals->setZero();
}
const std::vector<std::pair<double, Point>>& OneElectronIntegralController::getExternalCharges() {
  return _externalChargeController->getExternalCharges();
}
void OneElectronIntegralController::cacheExtChargeIntegrals(const MatrixInBasis<RESTRICTED>& integrals) {
  _extChargeIntegrals = std::make_unique<MatrixInBasis<RESTRICTED>>(integrals);
}
const MatrixInBasis<RESTRICTED>& OneElectronIntegralController::getExtChargeIntegrals() {
  if (!_extChargeIntegrals) {
    if (this->getExternalCharges().empty()) {
      throw SerenityError("Integrals over external charges were requested but no charges are available!");
    }
    calcExternalChargeIntegrals();
  }
  return *_extChargeIntegrals;
}
} /* namespace Serenity */
