/**
 * @file   Transformation.cpp
 *
 * @date   08-28-2013, 07-08-2016
 * @author Thomas Dresselhaus, Michael Boeckers
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
#include "basis/Transformation.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/wrappers/Libint.h"
#include "math/Matrix.h"
#include "math/linearAlgebra/Orthogonalization.h"
/* Include Std and External Headers */
#include <cassert>
#include <cmath>
#include <vector>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::unique_ptr<OrbitalController<SCFMode>> Transformation::transformMOs(OrbitalController<SCFMode>& orbitalsA,
                                                                         std::shared_ptr<BasisController> basisControllerB,
                                                                         const MatrixInBasis<RESTRICTED>& overlapB) {
  assert(basisControllerB);
  assert(overlapB.getBasisController() == basisControllerB);
  // coefficientsA of MOs in Basis A
  CoefficientMatrix<SCFMode> coefficientsA = orbitalsA.getCoefficients();
  // The basisController of Basis A
  auto basisControllerA = orbitalsA.getBasisController();
  // The number of AOs in Basis A
  const unsigned int nOrbitalsA = basisControllerA->getNBasisFunctions();
  // The number of AOs in Basis B
  const unsigned int nOrbitalsB = basisControllerB->getNBasisFunctions();

  // Mixed AO-overlap between basis A and basis B
  auto& libint = Libint::getInstance();
  auto overlapAB = libint.compute1eInts(LIBINT_OPERATOR::overlap, basisControllerA, basisControllerB);

  // Calculate inverse of AO overlap integrals in basis B. Use SVD,
  // i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd singularValues = svd.singularValues();
  for (unsigned int i = 0; i < singularValues.rows(); ++i) {
    // if singular value is too small, set it to zero
    if (singularValues(i) < 1.0e-10) {
      singularValues(i) = 0.0;
    }
    else {
      singularValues(i) = 1 / singularValues(i);
    }
  }
  // build inverse, i.e. V * D * U^T
  Eigen::MatrixXd overlapBinverse = svd.matrixV() * singularValues.asDiagonal() * svd.matrixU().transpose();
  // Calculate projection operator P_{BA}, i.e. the projector for the
  // transformation from basis A to basis B
  Eigen::MatrixXd projectionOperator = overlapBinverse * overlapAB;
  // Calculate new MO coefficients for basis B
  Eigen::MatrixXd newCoefficientsB = projectionOperator * coefficientsA;

  //  //Orthogonalize new coefficients
  //  Orthogonalization::modifiedGramSchmidt(newCoefficientsB);
  // Normalize new coefficients
  for (unsigned int i = 0; i < nOrbitalsA; i++) {
    double factor = 0;
    for (unsigned int j = 0; j < nOrbitalsB; j++) {
      for (unsigned int k = 0; k < nOrbitalsB; k++) {
        factor += newCoefficientsB(j, i) * newCoefficientsB(k, i) * overlapB(j, k);
      }
    }
    factor = sqrt(factor);
    for (unsigned int j = 0; j < nOrbitalsB; j++) {
      newCoefficientsB(j, i) /= factor;
    }
  }
  // Construct new orbitalSet
  auto coefficientsBptr = std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(basisControllerB));
  auto& coefficientsB = *coefficientsBptr;
  coefficientsB.setZero();
  coefficientsB.block(0, 0, nOrbitalsB, nOrbitalsA) = newCoefficientsB;
  // OrbitalController requires eigenvalues. These are set to zero.
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> eigenvaluesBptr(
      new SpinPolarizedData<SCFMode, Eigen::VectorXd>(nOrbitalsB));
  auto coreOrbitalPtr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(orbitalsA.getOrbitalFlags());
  // return new OrbitalController
  return std::make_unique<OrbitalController<SCFMode>>(std::move(coefficientsBptr), basisControllerB,
                                                      std::move(eigenvaluesBptr), std::move(coreOrbitalPtr));
}

template std::unique_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>>
Transformation::transformMOs(OrbitalController<Options::SCF_MODES::RESTRICTED>& orbitalsA,
                             std::shared_ptr<BasisController> basisControllerB, const MatrixInBasis<RESTRICTED>& overlapB);

} /* namespace Serenity */
