/**
 * @file   InitialGuessCalculator.cpp
 *
 * @date   Apr 22, 2015
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
#include "scf/initialGuess/InitialGuessCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cmath>

namespace Serenity {

void InitialGuessCalculator<Options::SCF_MODES::UNRESTRICTED>::scrambleOrbitals(
    OrbitalController<Options::SCF_MODES::UNRESTRICTED>& orbitals,
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nElectrons) {
  // Do not do anything for one-electron systems and other special cases
  if (nElectrons.alpha + nElectrons.beta <= 1 || nElectrons.alpha == 0 || nElectrons.beta == 0 ||
      nElectrons.alpha != nElectrons.beta)
    return;
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coefficients = orbitals.getCoefficients();
  auto& alphaCoefficients = coefficients.alpha;
  if (!(alphaCoefficients.cols() == coefficients.beta.cols()))
    throw SerenityError("Scramble orbitals: Different number of orbitals for alpha and beta.");
  const unsigned int nCol = alphaCoefficients.cols();
  const unsigned int nRow = alphaCoefficients.rows();
  const unsigned int nAlphaElec = nElectrons.alpha;
  // Cannot scramble with virtual orbitals if there are none
  if (nAlphaElec == nRow)
    return;
  const unsigned int scramble =
      (nRow < 4 || nAlphaElec <= 2 || nElectrons.beta <= 2 || nRow - nAlphaElec < 2 || nRow - nElectrons.beta < 2) ? 1 : 2;
  for (unsigned int i = 0; i < scramble; ++i) {
    for (unsigned int j = 0; j < nCol; ++j) {
      const auto sinus = sin(45.0 * M_PI / 180);
      const auto cosinus = cos(45.0 * M_PI / 180);
      const auto ij = alphaCoefficients(j, nAlphaElec - 1);
      const auto nColij = alphaCoefficients(j, nAlphaElec - 1 + i);
      alphaCoefficients(j, nAlphaElec - 1) = (cosinus * ij - sinus * nColij);
      alphaCoefficients(j, nAlphaElec - 1 + i) = cosinus * nColij + sinus * ij;
    }
    if (scramble == 2) {
      const auto sinus = sin(20.0 * M_PI / 180);
      const auto cosinus = cos(20.0 * M_PI / 180);
      for (unsigned int j = 0; j < nCol; ++j) {
        const auto ij = alphaCoefficients(j, nAlphaElec - 2);
        const auto nColij = alphaCoefficients(j, nAlphaElec + 1);
        alphaCoefficients(j, nAlphaElec - 2) = (cosinus * ij - sinus * nColij);
        alphaCoefficients(j, nAlphaElec + 1) = cosinus * nColij + sinus * ij;
      }
    }
  }
  orbitals.updateOrbitals(coefficients, orbitals.getEigenvalues());
}

std::unique_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>
UnrestrictedFromRestrictedGuess::calculateInitialGuess(std::shared_ptr<SystemController> systemController) {
  auto restrictedESguess = _restrictedGuessCalculator->calculateInitialGuess(systemController);
  auto restrictedOrbs = restrictedESguess->getMolecularOrbitals();
  // Create empty unrestricted orbital set
  auto unrestrictedOrbs = std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(
      new OrbitalController<Options::SCF_MODES::UNRESTRICTED>(restrictedOrbs->getBasisController(),
                                                              restrictedOrbs->getNCoreOrbitals()));
  // Copy orbitals to unrestricted
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> restrictedCoefficients = restrictedOrbs->getCoefficients();
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> unrestrictedCoefficients = unrestrictedOrbs->getCoefficients();
  unrestrictedCoefficients.alpha = restrictedCoefficients;
  unrestrictedCoefficients.beta = restrictedCoefficients;
  // Copy orbital energies to unrestricted
  SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd> restrictedEigenvalues = restrictedOrbs->getEigenvalues();
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd> unrestrictedEigenvalues =
      unrestrictedOrbs->getEigenvalues();
  unrestrictedEigenvalues.alpha = restrictedEigenvalues;
  unrestrictedEigenvalues.beta = restrictedEigenvalues;
  unrestrictedOrbs->updateOrbitals(unrestrictedCoefficients, unrestrictedEigenvalues);
  // Scramble orbitals
  this->scrambleOrbitals(*unrestrictedOrbs, systemController->getNElectrons<Options::SCF_MODES::UNRESTRICTED>());
  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> urES(
      new ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>(
          unrestrictedOrbs, restrictedESguess->getOneElectronIntegralController(),
          systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>()));

  return urES;
}

} // namespace Serenity
