/**
 * @file   HCoreGuessCalculator.cpp
 *
 * @date   Nov 7, 2013
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
#include "scf/initialGuess/HCoreGuessCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::unique_ptr<ElectronicStructure<SCFMode>>
HCoreGuessCalculator<SCFMode>::calculateInitialGuess(std::shared_ptr<SystemController> systemController) {
  /*
   * Get variables local
   */
  auto basisController = systemController->getBasisController();
  const unsigned int nOrbitals = basisController->getNBasisFunctions();
  auto settings = systemController->getSettings();

  /*
   * Create new set of orbitals
   */
  auto eigenvalues = std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
      new SpinPolarizedData<SCFMode, Eigen::VectorXd>(nOrbitals));
  auto& eps = *eigenvalues;
  auto coefficientMatrix = std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(basisController));
  auto& c = *coefficientMatrix;
  const auto oneIntController = systemController->getOneElectronIntegralController();
  /*
   * Solve h*C= S*C*eps
   */
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(oneIntController->getOneElectronIntegrals(),
                                                               oneIntController->getOverlapIntegrals());
  for_spin(c) {
    c_spin = es.eigenvectors();
  };
  for_spin(eps) {
    for (unsigned int i = 0; i < nOrbitals; ++i) {
      eps_spin[i] = es.eigenvalues()[i];
    }
  };

  auto orbs =
      std::make_shared<OrbitalController<SCFMode>>(std::move(coefficientMatrix), systemController->getBasisController(),
                                                   *eigenvalues, systemController->getNCoreElectrons());
  orbs->setCanOrthTh(systemController->getSettings().scf.canOrthThreshold);
  std::unique_ptr<ElectronicStructure<SCFMode>> elecStruct(new ElectronicStructure<SCFMode>(
      orbs, systemController->getOneElectronIntegralController(), systemController->getNOccupiedOrbitals<SCFMode>()));

  // Use fractional occupations if requested.
  elecStruct->getDensityMatrixController()->setDegeneracyThreshold(settings.scf.degeneracyThreshold);

  return elecStruct;
}

template class HCoreGuessCalculator<Options::SCF_MODES::RESTRICTED>;
template class HCoreGuessCalculator<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
