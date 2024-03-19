/**
 * @file ROHF.cpp
 *
 * @date Jun 21, 2023
 * @author Niklas Niemeyer
 *
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
#include "scf/ROHF.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/OneElectronIntegralController.h"
#include "misc/SerenityError.h"

namespace Serenity {

template<>
void ROHF<Options::SCF_MODES::RESTRICTED>::addConstraint(FockMatrix<RESTRICTED>& F,
                                                         std::shared_ptr<ElectronicStructure<RESTRICTED>> es,
                                                         Options::ROHF_TYPES rohf, double suhfLambda) {
  (void)F;
  (void)es;
  (void)rohf;
  (void)suhfLambda;
  return;
}

template<>
void ROHF<Options::SCF_MODES::UNRESTRICTED>::addConstraint(FockMatrix<UNRESTRICTED>& F,
                                                           std::shared_ptr<ElectronicStructure<UNRESTRICTED>> es,
                                                           Options::ROHF_TYPES rohf, double suhfLambda) {
  const auto& S = es->getOneElectronIntegralController()->getOverlapIntegrals();
  DensityMatrix<UNRESTRICTED> dmat = es->getDensityMatrix();
  if (rohf == Options::ROHF_TYPES::CUHF) {
    Eigen::MatrixXd alphaCoeff = es->getMolecularOrbitals()->getCoefficients().alpha;
    Eigen::MatrixXd alphaCoeffinv = alphaCoeff.transpose() * S;

    // Choose an orthogonal basis to represent the orbitals in (alpha MOs).
    Eigen::MatrixXd totalDensity = alphaCoeffinv * dmat.total() * alphaCoeffinv.transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(totalDensity);
    Eigen::MatrixXd noCoeff = alphaCoeff * eigensolver.eigenvectors();
    Eigen::MatrixXd deltaUHF = -0.5 * noCoeff.transpose() * F.difference() * noCoeff;

    unsigned noa = es->getNOccupiedOrbitals().alpha;
    unsigned nob = es->getNOccupiedOrbitals().beta;

    unsigned nCore = std::min(noa, nob);
    unsigned nActive = std::max(noa, nob) - nCore;
    unsigned nVirtual = S.rows() - nActive - nCore;

    // Null all blocks but the core-virtual and virtual-core one.
    deltaUHF.block(0, 0, nVirtual + nActive, nVirtual + nActive).setZero();
    deltaUHF.block(nVirtual, nVirtual, nActive + nCore, nActive + nCore).setZero();

    Eigen::MatrixXd noCoeffinv = noCoeff.transpose() * S;

    deltaUHF = (noCoeffinv.transpose() * deltaUHF * noCoeffinv).eval();

    F.alpha += deltaUHF;
    F.beta -= deltaUHF;
  }
  else if (rohf == Options::ROHF_TYPES::SUHF) {
    F.alpha -= 2.0 * suhfLambda * S * dmat.beta * S;
    F.beta -= 2.0 * suhfLambda * S * dmat.alpha * S;
  }
  else {
    throw SerenityError("ROHF type not supported");
  }
}

template class ROHF<Options::SCF_MODES::RESTRICTED>;
template class ROHF<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
