/**
 * @file SCFAnalysis.cpp
 *
 * @date Dec 1, 2016
 * @author M. Boeckers
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
#include "scf/SCFAnalysis.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/MatrixInBasis.h"
#include "energies/EnergyContributions.h"
#include "integrals/wrappers/Libint.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES T>
SCFAnalysis<T>::SCFAnalysis(std::shared_ptr<SystemController> systemController,
                            std::shared_ptr<OneElectronIntegralController> oneIntController,
                            std::shared_ptr<EnergyComponentController> energyController)
  : _systemController(systemController), _oneIntController(oneIntController), _energyController(energyController) {
}

template<>
double SCFAnalysis<Options::SCF_MODES::RESTRICTED>::S2() {
  return 0.0;
}

template<>
double SCFAnalysis<Options::SCF_MODES::UNRESTRICTED>::S2() {
  double S = 0.5 * _systemController->getSpin();
  if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    //
    // Calculate <S*S> as functional of density
    //

    // Calculate spin density
    auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        _systemController->getSettings(), _systemController->getBasisController(), _systemController->getGridController());
    auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
        basisFunctionOnGridController, _systemController->getSettings().grid.blockAveThreshold);
    auto densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
        _systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>(),
        _systemController->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
    auto densityOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(
        densityOnGridCalculator, densityMatrixController);
    auto& densityOnGrid = densityOnGridController->getDensityOnGrid();
    Eigen::VectorXd spinDensity = densityOnGrid.alpha - densityOnGrid.beta;

    // Integrage over all regions where spin density is negative
    for (unsigned int i = 0; i < spinDensity.rows(); ++i) {
      if (spinDensity(i) > 0)
        spinDensity(i) = 0.0;
    }
    double integral = spinDensity.dot(_systemController->getGridController()->getWeights());

    return S * (S + 1) - integral;
  }
  else if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    //
    // Calculate <S*S> from UHF orbitals
    //
    auto overlap = _oneIntController->getOverlapIntegrals();
    auto nEl = _systemController->getNElectrons<Options::SCF_MODES::UNRESTRICTED>();
    auto coeff = _systemController->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
    Eigen::MatrixXd Sab = coeff.alpha.block(0, 0, overlap.rows(), nEl.alpha).transpose() * overlap *
                          coeff.beta.block(0, 0, overlap.rows(), nEl.beta);
    return S * (S + 1) + nEl.beta - Sab.array().square().sum();
  }
  throw SerenityError("SCFAnalysis: No analysis available for the requested method.");
}

template<Options::SCF_MODES T>
double SCFAnalysis<T>::VirialRatio() {
  // Calculate kinetic energy
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd TInts = libint.compute1eInts(libint2::Operator::kinetic, _systemController->getBasisController());
  auto P = _systemController->getElectronicStructure<T>()->getDensityMatrix();
  double TEnergy = 0.0;
  for_spin(P) {
    TEnergy += TInts.cwiseProduct(P_spin).sum();
  };
  // Calculate potential energy
  double VEnergy = 0.0;
  if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    VEnergy = _energyController->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY) - TEnergy;
  }
  else if (_systemController->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    VEnergy = _energyController->getEnergyComponent(ENERGY_CONTRIBUTIONS::HF_ENERGY) - TEnergy;
  }
  else {
    assert(false);
  }

  return -1.0 * VEnergy / TEnergy;
}

template class SCFAnalysis<Options::SCF_MODES::RESTRICTED>;
template class SCFAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
