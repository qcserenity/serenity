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
#include "data/grid/DensityOnGridFactory.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "data/matrices/MatrixInBasis.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/wrappers/Libint.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SCFAnalysis<SCFMode>::SCFAnalysis(std::vector<std::shared_ptr<SystemController>> systemController,
                                  std::shared_ptr<GridController> supersystemGrid)
  : _systemController(systemController), _supersystemGrid(supersystemGrid) {
}

template<>
double SCFAnalysis<Options::SCF_MODES::RESTRICTED>::getS2(bool useUHForbitals) {
  (void)useUHForbitals; // no warnings
  return 0.0;
}

template<>
double SCFAnalysis<Options::SCF_MODES::UNRESTRICTED>::getS2(bool useUHForbitals) {
  /*
   * Calculate physical <S*S>
   */
  double S = 0.0;
  double integral = 0.0;
  if (_systemController.size() == 1) {
    S = 0.5 * _systemController[0]->getSpin();
  }
  else {
    for (unsigned int i = 0; i < _systemController.size(); i++) {
      S += 0.5 * _systemController[i]->getSpin();
    }
  }

  /*
   * Calculate <S*S> as functional of density
   */
  if (_systemController[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT and useUHForbitals == false) {
    // Calculate spin density
    std::shared_ptr<DensityOnGridController<Options::SCF_MODES::UNRESTRICTED>> densityOnGridController(nullptr);
    /*
     * Supersystem Calculation
     */
    if (_systemController.size() == 1) {
      auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
          _systemController[0]->getSettings(), _systemController[0]->getBasisController(),
          _systemController[0]->getGridController());
      auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
          basisFunctionOnGridController, _systemController[0]->getSettings().grid.blockAveThreshold);
      auto densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
          _systemController[0]->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>(),
          _systemController[0]->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
      densityOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(
          densityOnGridCalculator, densityMatrixController);
      auto& densityOnGrid = densityOnGridController->getDensityOnGrid();
      Eigen::VectorXd spinDensity = densityOnGrid.alpha - densityOnGrid.beta;

      // Integrate over all regions where spin density is negative
      for (unsigned int i = 0; i < spinDensity.rows(); ++i) {
        if (spinDensity(i) > 0)
          spinDensity(i) = 0.0;
      }
      integral = spinDensity.dot(_systemController[0]->getGridController()->getWeights());
    }
    /*
     * Subsystem Calculation
     */
    else {
      if (_supersystemGrid == nullptr) {
        auto superSystemGeometry = std::make_shared<Geometry>();
        for (unsigned int i = 0; i < _systemController.size(); i++) {
          *superSystemGeometry += *_systemController[i]->getGeometry();
        }
        superSystemGeometry->deleteIdenticalAtoms();
        // supersystem grid
        Options::GRID_PURPOSES gridacc = Options::GRID_PURPOSES::DEFAULT;
        _supersystemGrid = AtomCenteredGridControllerFactory::produce(superSystemGeometry,
                                                                      _systemController[0]->getSettings().grid, gridacc);
      }

      std::vector<std::shared_ptr<DensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>> _densOnGridControllers;
      for (unsigned int i = 0; i < _systemController.size(); i++) {
        auto dens = DensityOnGridFactory<Options::SCF_MODES::UNRESTRICTED>::produce(
            _systemController[i]->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrixController(),
            _supersystemGrid, 1, _systemController[0]->getSettings());
        _densOnGridControllers.push_back(dens);
      }
      densityOnGridController =
          std::make_shared<SupersystemDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(_densOnGridControllers);
      auto& densityOnGrid = densityOnGridController->getDensityOnGrid();
      Eigen::VectorXd spinDensity = densityOnGrid.alpha - densityOnGrid.beta;

      // Integrate over all regions where spin density is negative
      for (unsigned int i = 0; i < spinDensity.rows(); ++i) {
        if (spinDensity(i) > 0)
          spinDensity(i) = 0.0;
      }
      integral = spinDensity.dot(_supersystemGrid->getWeights());
    }

    return S * (S + 1) - integral;

    /*
     * Calculate <S*S> from UHF orbitals
     */
  }
  else if (_systemController[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF or useUHForbitals == true) {
    /*
     * Supersystem Calculation
     */
    if (_systemController.size() == 1) {
      auto overlap = _systemController[0]->getOneElectronIntegralController()->getOverlapIntegrals();
      auto nEl = _systemController[0]->getNElectrons<Options::SCF_MODES::UNRESTRICTED>();
      auto coeff = _systemController[0]->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
      Eigen::MatrixXd Sab = coeff.alpha.block(0, 0, overlap.rows(), nEl.alpha).transpose() * overlap *
                            coeff.beta.block(0, 0, overlap.rows(), nEl.beta);
      return S * (S + 1) + nEl.beta - Sab.array().square().sum();
    }
    else {
      throw SerenityError(
          "No method implemented for the calculation of <S*S> from UHF orbitals for a subsystem approach");
    }
  }
  else {
    throw SerenityError("SCFAnalysis: No analysis available for the requested method.");
  }
}

template<Options::SCF_MODES SCFMode>
double SCFAnalysis<SCFMode>::getVirialRatio() {
  if (_systemController.size() == 1) {
    // Calculate kinetic energy
    auto& libint = Libint::getInstance();
    Eigen::MatrixXd TInts = libint.compute1eInts(LIBINT_OPERATOR::kinetic, _systemController[0]->getBasisController());
    auto P = _systemController[0]->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    double TEnergy = 0.0;
    for_spin(P) {
      TEnergy += TInts.cwiseProduct(P_spin).sum();
    };
    // Calculate potential energy
    double VEnergy = 0.0;
    if (_systemController[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
      VEnergy = _systemController[0]->template getElectronicStructure<SCFMode>()->getEnergyComponentController()->getEnergyComponent(
                    ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY) -
                TEnergy;
    }
    else if (_systemController[0]->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
      VEnergy = _systemController[0]->template getElectronicStructure<SCFMode>()->getEnergyComponentController()->getEnergyComponent(
                    ENERGY_CONTRIBUTIONS::HF_ENERGY) -
                TEnergy;
    }
    else {
      throw SerenityError("Unknown electronic structure method.");
    }
    return -1.0 * VEnergy / TEnergy;
  }
  else {
    throw SerenityError("Virial ratio not implemented for more than one system.");
    return 0.0;
  }
}

template class SCFAnalysis<Options::SCF_MODES::RESTRICTED>;
template class SCFAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
