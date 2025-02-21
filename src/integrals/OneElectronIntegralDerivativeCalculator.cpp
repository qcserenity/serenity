/**
 * @file   OneElectronIntegralDerivativeCalculator.cpp
 *
 * @date   Jun 17, 2024
 * @author Anton Rikus
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "integrals/OneElectronIntegralDerivativeCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h"                   // Loop shells
#include "data/ExternalChargeController.h" // External point charges.
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Geometry.h"
#include "geometry/Point.h"
#include "integrals/wrappers/Libint.h"
#include "math/linearAlgebra/MatrixFunctions.h" // symmetrize
#include "misc/Timing.h"                        // Timings

namespace Serenity {

OneElectronIntegralDerivativeCalculator::OneElectronIntegralDerivativeCalculator(
    std::shared_ptr<AtomCenteredBasisController> basisController, std::shared_ptr<const Geometry> geometry,
    const std::vector<std::pair<double, Point>> pointCharges)
  : _basisController(basisController), _geometry(geometry), _pointCharges(pointCharges) {
}

OneElectronIntegralDerivativeCalculator::OneElectronIntegralDerivativeCalculator(
    std::shared_ptr<AtomCenteredBasisController> basisController, std::shared_ptr<const Geometry> geometry,
    std::shared_ptr<ExternalChargeController> externalChargeController)
  : OneElectronIntegralDerivativeCalculator(basisController, geometry, getAllCharges(geometry, externalChargeController)) {
}

const Eigen::MatrixXd OneElectronIntegralDerivativeCalculator::getNucKinDerivative(const MatrixInBasis<RESTRICTED>& density) {
  Timings::takeTime("Derivatives -   1e-Integrals");
  Eigen::MatrixXd symDens = symmetrize(density);

  unsigned nAtoms = _pointCharges.empty() ? this->_geometry->getNAtoms() : _pointCharges.size();
  auto mapping = _basisController->getAtomIndicesOfBasisShells();
  auto& basis = _basisController->getBasis();

  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);
  Libint& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::kinetic, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(this->_geometry->getNAtoms(), 3);
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < basis.size(); i++) {
      unsigned int offI = _basisController->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; j++) {
        unsigned int offJ = _basisController->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        if (libint.compute(LIBINT_OPERATOR::kinetic, 1, *basis[i], *basis[j], intDerivs)) {
          for (unsigned int iCol = 0; iCol < 2; ++iCol) {
            unsigned iAtom = iCol ? mapping[j] : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) +=
                  ((i == j) ? 1.0 : 2.0) * tmp.cwiseProduct(symDens.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }
#pragma omp critical
    { gradientContr.topRows(this->_geometry->getNAtoms()) += gradientContrPriv; }
  } /* END OpenMP parallel */

  if (_pointCharges.empty()) {
    libint.initialize(LIBINT_OPERATOR::nuclear, 1, 2, _geometry->getAtoms());
  }
  else {
    libint.initialize(LIBINT_OPERATOR::nuclear, 1, 2, _pointCharges);
  }

#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < basis.size(); i++) {
      unsigned int offI = _basisController->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; j++) {
        unsigned int offJ = _basisController->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        if (libint.compute(LIBINT_OPERATOR::nuclear, 1, *basis[i], *basis[j], intDerivs)) {
          for (unsigned int iCol = 0; iCol < (2 + nAtoms); ++iCol) {
            unsigned iAtom = iCol ? ((iCol > 1) ? iCol - 2 : mapping[j]) : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) +=
                  ((i == j) ? 1.0 : 2.0) * tmp.cwiseProduct(symDens.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }
#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  libint.finalize(LIBINT_OPERATOR::kinetic, 1, 2);
  libint.finalize(LIBINT_OPERATOR::nuclear, 1, 2);

  Timings::timeTaken("Derivatives -   1e-Integrals");
  return gradientContr;
}

const Eigen::MatrixXd OneElectronIntegralDerivativeCalculator::getOverlapDerivative(const MatrixInBasis<RESTRICTED>& density) {
  Timings::takeTime("Derivatives -   1e-Integrals");
  Eigen::MatrixXd symDens = symmetrize(density);

  unsigned nAtoms = this->_geometry->getNAtoms();
  auto mapping = _basisController->getAtomIndicesOfBasisShells();
  auto& basis = _basisController->getBasis();

  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);
  Libint& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::overlap, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < basis.size(); i++) {
      unsigned int offI = _basisController->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; j++) {
        unsigned int offJ = _basisController->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        if (libint.compute(LIBINT_OPERATOR::overlap, 1, *basis[i], *basis[j], intDerivs)) {
          for (unsigned int iCol = 0; iCol < 2; ++iCol) {
            unsigned iAtom = iCol ? mapping[j] : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) +=
                  ((i == j) ? 1.0 : 2.0) * tmp.cwiseProduct(symDens.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }
#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  libint.finalize(LIBINT_OPERATOR::overlap, 1, 2);

  Timings::timeTaken("Derivatives -   1e-Integrals");
  return gradientContr;
}

std::vector<std::pair<double, Point>>
OneElectronIntegralDerivativeCalculator::getAllCharges(std::shared_ptr<const Geometry> geometry,
                                                       std::shared_ptr<ExternalChargeController> externalChargeController) {
  std::vector<std::pair<double, Point>> allCharges;
  for (const auto& atom : geometry->getAtoms()) {
    allCharges.emplace_back(atom->getEffectiveCharge(), *atom);
  }
  const auto extCharges = externalChargeController->getExternalCharges();
  allCharges.insert(allCharges.end(), extCharges.begin(), extCharges.end());
  return allCharges;
}

} /* namespace Serenity */