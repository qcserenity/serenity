/**
 * @file SPADEAlgorithm.cpp
 *
 * @date 8 Mar 2020
 * @author Moritz Bensberg
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
#include "analysis/orbitalLocalization/SPADEAlgorithm.h"
/* Include Serenity Internal Headers */
#include "basis/BasisFunctionMapper.h"               //Block selection in coefficient matrix.
#include "data/OrbitalController.h"                  //Coefficients and orbital update.
#include "data/matrices/CoefficientMatrix.h"         //Coefficient matrix.
#include "data/matrices/MatrixInBasis.h"             //Overlap matrix.
#include "integrals/OneElectronIntegralController.h" //Overlap integrals for symmetric orthogonalization.
#include "io/FormattedOutputStream.h"                //Filtered output streams.
#include "math/linearAlgebra/MatrixFunctions.h"      //Matrix sqrt
#include "misc/SerenityError.h"                      //Error messages
#include "system/SystemController.h"                 //System controller defintion.
namespace Serenity {

template<Options::SCF_MODES SCFMode>
SPADEAlgorithm<SCFMode>::SPADEAlgorithm(std::shared_ptr<SystemController> supersystem, std::shared_ptr<SystemController> activeSystem)
  : _supersystem(supersystem), _activeSystem(activeSystem) {
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXi> SPADEAlgorithm<SCFMode>::run() {
  auto orbitalController = _supersystem->getActiveOrbitalController<SCFMode>();
  const MatrixInBasis<RESTRICTED>& S = _supersystem->getOneElectronIntegralController()->getOverlapIntegrals();
  const Eigen::MatrixXd sqrtS = mSqrt_Sym(S);
  CoefficientMatrix<SCFMode> coefficients = orbitalController->getCoefficients();
  auto nOcc = _supersystem->getNOccupiedOrbitals<SCFMode>();
  SpinPolarizedData<SCFMode, Eigen::VectorXi> assignment;
  BasisFunctionMapper basisFunctionMapper(_supersystem->getBasisController());
  auto projection = basisFunctionMapper.getSparseProjection(_activeSystem->getBasisController());
  for_spin(coefficients, nOcc, assignment) {
    const Eigen::MatrixXd oldC = coefficients_spin.leftCols(nOcc_spin).eval();
    const Eigen::MatrixXd d = oldC * oldC.transpose();
    const Eigen::MatrixXd orthoC = sqrtS * oldC;
    const Eigen::MatrixXd orthoC_A = *projection * orthoC;
    // Calculation of the full matrix V is necessary for the orbital rotation.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(orthoC_A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::MatrixXd rightSingularVectors = svd.matrixV();
    const Eigen::VectorXd singularValues = svd.singularValues();
    const Eigen::MatrixXd newC = (oldC * rightSingularVectors).eval();
    const Eigen::MatrixXd d2 = newC * newC.transpose();
    if ((d - d2).array().abs().sum() > 1e-10)
      throw SerenityError("Density matrix changed during SPADE orbital construction.");
    coefficients_spin.leftCols(nOcc_spin) = newC;
    unsigned int lastActiveOrbital = nOcc_spin;
    double largestDifference = 0.0;
    for (unsigned int iOrb = 1; iOrb < singularValues.size(); iOrb++) {
      double singularValueDifference = std::fabs(singularValues[iOrb - 1] - singularValues[iOrb]);
      if (singularValueDifference > largestDifference) {
        largestDifference = singularValueDifference;
        lastActiveOrbital = iOrb;
      }
    } // for iOrb
    assignment_spin = Eigen::VectorXi::Constant(nOcc_spin, 1);
    assignment_spin.head(lastActiveOrbital) = Eigen::VectorXi::Zero(lastActiveOrbital);
  };
  // TODO: Allow the splitting of core and non-core orbitals.
  orbitalController->updateOrbitals(coefficients, _supersystem->getActiveOrbitalController<SCFMode>()->getEigenvalues());
  return assignment;
}

template class SPADEAlgorithm<Options::SCF_MODES::RESTRICTED>;
template class SPADEAlgorithm<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
