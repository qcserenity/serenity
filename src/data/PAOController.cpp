/**
 * @file PAOController.cpp
 *
 * @date Jan 27, 2019
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
#include "data/PAOController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisFunctionMapper.h"
#include "io/FormattedOutputStream.h" //Filtered output streams.

namespace Serenity {
void PAOController::constructPAOs() {
  if (_paoCoefficients)
    return;
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto activeBasisController = _activeDensity->getBasisController();
  _paoCoefficients.reset(new CoefficientMatrix<scfMode>(activeBasisController));
  unsigned int nBasisFunctions = activeBasisController->getNBasisFunctions();
  auto& coefficients = *_paoCoefficients;
  coefficients = Eigen::MatrixXd::Identity(nBasisFunctions, nBasisFunctions);
  // Active system
  const auto& d = *_activeDensity;
  const auto& s = *_s;
  // Environment - removing projected occupied orbitals
  Eigen::MatrixXd d_env = Eigen::MatrixXd::Zero(nBasisFunctions, nBasisFunctions);
  BasisFunctionMapper mapper(_activeDensity->getBasisController());
  for (auto& envDensity : _environmentDensities) {
    auto projection = mapper.getSparseProjection(envDensity->getBasisController());
    d_env += projection->transpose() * (*envDensity) * (*projection);
  }
  coefficients -= 0.5 * (d + d_env) * s;
  // Search for functions which are no longer present in the basis.
  // Otherwise renormalize.
  Eigen::MatrixXd metric = coefficients.transpose() * s * coefficients;
  for (unsigned int iCol = 0; iCol < nBasisFunctions; ++iCol) {
    double norm = sqrt(metric(iCol, iCol));
    if (norm < _paoNormalizationThreshold) {
      OutputControl::dOut << "Vanishing PAO for basis function " << iCol << std::endl;
      coefficients.col(iCol).setZero();
    }
    else {
      coefficients.col(iCol) *= 1.0 / norm;
    }
  }
}

const MatrixInBasis<RESTRICTED>& PAOController::getS_PAO() {
  if (!_s_pao) {
    _s_pao = std::make_unique<MatrixInBasis<RESTRICTED>>(_s->getBasisController());
    const auto& paoCoefficients = this->getAllPAOs();
    *_s_pao = paoCoefficients.transpose() * *_s * paoCoefficients;
  }
  return *_s_pao;
}

} /* namespace Serenity */
