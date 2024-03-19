/**
 * @file   ExtendedHueckel.cpp
 *
 * @date   Nov 23, 2013
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
#include "scf/initialGuess/ExtendedHueckel.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h"
#include "basis/Transformation.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

const std::map<std::string, std::vector<double>> ExtendedHueckel::_parameters = ExtendedHueckel::assignParameters();

std::unique_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>>
ExtendedHueckel::calculateHueckelOrbitals(const std::shared_ptr<SystemController> systemController,
                                          std::shared_ptr<AtomCenteredBasisController> minimalBasisController,
                                          const MatrixInBasis<RESTRICTED>& minimalBasisOverlaps) {
  // Abbreviations
  const auto& geometry = systemController->getGeometry();
  const auto& atoms = geometry->getAtoms();
  const auto& minimalBasis = minimalBasisController->getBasis();
  const auto& atomBasisIndices = minimalBasisController->getBasisIndicesRed();
  const unsigned int nBasisFunctions = minimalBasisController->getNBasisFunctions();
  /*
   * Construct Hückel matrix
   */
  MatrixInBasis<RESTRICTED> hueckelMatrix(minimalBasisController);
  /*
   * Fill diagonal elements with tabulated values.
   */
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    if (!((_parameters.at(atoms[i]->getAtomType()->getElementSymbol())).size() ==
          atoms[i]->getNBasisFunctions(minimalBasisController->getBasisLabel())))
      throw SerenityError("No hueckel parameters known for atom '" + atoms[i]->getAtomType()->getElementSymbol() + "'.");
    const auto& paramsForAtom = _parameters.at(atoms[i]->getAtomType()->getElementSymbol());
    // Loop over basis function shells of this atom; j is an index for the basis function vector
    for (unsigned int j = atomBasisIndices[i].first; j != atomBasisIndices[i].second; ++j) {
      /*
       * In this loop cycle we fill the matrix starting at the first index belonging to this shell
       * of basis functions
       */
      const auto& firstIndexShell = minimalBasisController->extendedIndex(j);
      /*
       * The atom parameter vector of course only has entries for basis functions of that atom, thus
       * we need to subtract the first index we use for accessing the basis function vector.
       */
      const auto& thisParam = paramsForAtom[j - atomBasisIndices[i].first];
      // Loop over the basis functions within a shell
      for (unsigned int k = 0; k < minimalBasis[j]->getNContracted(); ++k) {
        hueckelMatrix(firstIndexShell + k, firstIndexShell + k) = thisParam;
      }
    }
  }
  /*
   * Calculate off-diagonal elements.
   */
  for (unsigned int i = 1; i < nBasisFunctions; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      hueckelMatrix(i, j) = 0.5 * _k * (hueckelMatrix(i, i) + hueckelMatrix(j, j)) * minimalBasisOverlaps(i, j);
    }
  }
  /*
   * Solve the eigenvalue equation
   */
  std::unique_ptr<SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>> eigenvalues(
      new SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>(nBasisFunctions));
  std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::RESTRICTED>> coefficients(
      new CoefficientMatrix<Options::SCF_MODES::RESTRICTED>(minimalBasisController));

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(hueckelMatrix, minimalBasisOverlaps);
  auto& c = *coefficients;
  for_spin(c) {
    c_spin = es.eigenvectors();
  };
  for (unsigned int i = 0; i != eigenvalues->size(); ++i) {
    (*eigenvalues)[i] = es.eigenvalues()[i];
  }
  /*
   * Construct and return set of final orbitals
   */
  return std::make_unique<OrbitalController<Options::SCF_MODES::RESTRICTED>>(
      std::move(coefficients), minimalBasisController, *eigenvalues, systemController->getNCoreElectrons());
}

std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
ExtendedHueckel::calculateInitialGuess(const std::shared_ptr<SystemController> systemController) {
  /*
   * Collect necessary data from the systemController and create the Hueckel orbitals
   */
  auto minimalBasisController = systemController->getAtomCenteredBasisController(Options::BASIS_PURPOSES::HUECKEL);
  FockMatrix<Options::SCF_MODES::RESTRICTED> minimalBasisOverlaps(minimalBasisController);
  minimalBasisOverlaps =
      systemController->getOneElectronIntegralController(Options::BASIS_PURPOSES::HUECKEL)->getOverlapIntegrals();
  auto hueckelOrbs = calculateHueckelOrbitals(systemController, minimalBasisController, minimalBasisOverlaps);
  /*
   * Transform to the basis which is actually used
   */
  std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>> guessOrbs(
      Transformation::transformMOs(*hueckelOrbs, systemController->getBasisController(),
                                   systemController->getOneElectronIntegralController()->getOverlapIntegrals()));
  guessOrbs->setCanOrthTh(systemController->getSettings().scf.canOrthThreshold);
  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> elecStruct(
      new ElectronicStructure<Options::SCF_MODES::RESTRICTED>(
          guessOrbs, systemController->getOneElectronIntegralController(),
          systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>()));

  return elecStruct;
}

} /* namespace Serenity */
