/**
 * @file   RIIntegralDerivativeCalculator.cpp
 *
 * @date   Jan 10, 2025
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
#include "integrals/RIIntegralDerivativeCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h" // Loop shells
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Geometry.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "math/linearAlgebra/MatrixFunctions.h" // symmetrize
#include "misc/Timing.h"                        // Timings

namespace Serenity {

RIIntegralDerivativeCalculator::RIIntegralDerivativeCalculator(std::shared_ptr<AtomCenteredBasisController> basisController,
                                                               std::shared_ptr<AtomCenteredBasisController> auxBasisController)
  : _basisController(basisController), _auxBasisController(auxBasisController) {
}

const Eigen::MatrixXd
RIIntegralDerivativeCalculator::getCoulombIntegralDerivatives(const MatrixInBasis<RESTRICTED>& density1,
                                                              const MatrixInBasis<RESTRICTED>& density2a,
                                                              const MatrixInBasis<RESTRICTED>& density2b) {
  // X
  Eigen::MatrixXd d1 = symmetrize(density1);
  // D
  Eigen::MatrixXd d2a = symmetrize(density2a);
  // P (all in the case of TDDFT)
  Eigen::MatrixXd d2b = symmetrize(density2b);

  unsigned nAtoms = this->_basisController->getGeometry().getNAtoms();
  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);

  unsigned nXBasis = _auxBasisController->getNBasisFunctions();
  Eigen::MatrixXd sumMat1 = Eigen::MatrixXd::Zero(nXBasis, omp_get_max_threads());
  Eigen::MatrixXd sumMat2a = Eigen::MatrixXd::Zero(nXBasis, omp_get_max_threads());
  Eigen::MatrixXd sumMat2b = Eigen::MatrixXd::Zero(nXBasis, omp_get_max_threads());

  auto distribute1 = [&](unsigned p, unsigned q, unsigned K, double integral, unsigned threadId) {
    double perm = (p == q ? 1.0 : 2.0);
    sumMat1(K, threadId) += perm * integral * d1(p, q);
    sumMat2a(K, threadId) += perm * integral * d2a(p, q);
    sumMat2b(K, threadId) += perm * integral * d2b(p, q);
  };

  std::shared_ptr<RI_J_IntegralController> riInts =
      RI_J_IntegralControllerFactory::getInstance().produce(_basisController, _auxBasisController);
  riInts->loopOver3CInts(distribute1);

  Eigen::VectorXd coeff1 = riInts->getLLTMetric().solve(sumMat1.rowwise().sum()).eval();
  Eigen::VectorXd coeff2a = riInts->getLLTMetric().solve(sumMat2a.rowwise().sum()).eval();
  Eigen::VectorXd coeff2b = riInts->getLLTMetric().solve(sumMat2b.rowwise().sum()).eval();

  std::vector<unsigned int> atomIndicesOfBasis = _basisController->getAtomIndicesOfBasis();
  std::vector<unsigned int> atomIndicesOfAuxBasis = _auxBasisController->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
  auto const distribute2 = [&](const unsigned int& p, const unsigned int& q, const unsigned int& K,
                               Eigen::VectorXd& intValues, const unsigned int threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfAuxBasis[K], atomIndicesOfBasis[p], atomIndicesOfBasis[q]};
    double perm = (p == q ? 1.0 : 2.0);
    for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
      for (unsigned int iAtom = 0; iAtom < 3; ++iAtom) {
        coulombGradientPriv[threadID](atomMaps[iAtom], iDirection) +=
            perm * intValues(iAtom * 3 + iDirection) *
            (2 * coeff1(K) * d1(p, q) + coeff2b(K) * d2a(p, q) + coeff2a(K) * d2b(p, q));
      }
    }
  };
  TwoElecThreeCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _basisController, _auxBasisController,
                                     _basisController->getPrescreeningThreshold());
  looper.loop(distribute2);
  for (Eigen::MatrixXd& deriv : coulombGradientPriv) {
    gradientContr += deriv;
    deriv.setZero();
  }

  const auto& auxBasis = _auxBasisController->getBasis();
  const auto auxShellMapping = _auxBasisController->getAtomIndicesOfBasisShells();

  Libint& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::coulomb, 1, 2);
#pragma omp parallel for schedule(static, 1)
  for (unsigned int auxJ = 0; auxJ < auxBasis.size(); auxJ++) {
    const unsigned int nAuxJ = auxBasis[auxJ]->getNContracted();
    unsigned int offAuxJ = _auxBasisController->extendedIndex(auxJ);
    for (unsigned int auxK = 0; auxK <= auxJ; auxK++) {
      const unsigned int threadId = omp_get_thread_num();
      Eigen::MatrixXd intDerivs;
      const unsigned int nAuxK = auxBasis[auxK]->getNContracted();
      unsigned int offAuxK = _auxBasisController->extendedIndex(auxK);

      if (libint.compute(LIBINT_OPERATOR::coulomb, 1, *auxBasis[auxJ], *auxBasis[auxK], intDerivs)) {
        // prefac has nAuxJ rows and nAuxK columns
        Eigen::MatrixXd prefac = coeff2a.segment(offAuxJ, nAuxJ) * coeff2b.segment(offAuxK, nAuxK).transpose() +
                                 coeff1.segment(offAuxJ, nAuxJ) * coeff1.segment(offAuxK, nAuxK).transpose();
        prefac *= (auxJ == auxK ? 1.0 : 2.0);

        std::vector<unsigned> atomMaps = {auxShellMapping[auxJ], auxShellMapping[auxK]};
        for (unsigned int iAtom = 0; iAtom < 2; iAtom++) {
          for (unsigned int iDirection = 0; iDirection < 3; iDirection++) {
            Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nAuxK, nAuxJ);
            coulombGradientPriv[threadId](atomMaps[iAtom], iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
          }
        }
      }
    }
  }
  libint.finalize(LIBINT_OPERATOR::coulomb, 1, 2);

  for (const Eigen::MatrixXd& deriv : coulombGradientPriv) {
    gradientContr += deriv;
  }

  return gradientContr;
}

} /* namespace Serenity */