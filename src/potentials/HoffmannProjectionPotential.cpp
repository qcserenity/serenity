/**
 * @file HoffmannProjectionPotential.cpp
 *
 * @author Moritz Bensberg
 * @date Aug 26, 2019
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
#include "potentials/HoffmannProjectionPotential.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/OneElectronIntegralController.h"
#include "misc/SystemSplittingTools.h"
#include "potentials/HuzinagaProjectionPotential.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HoffmannProjectionPotential<SCFMode>::HoffmannProjectionPotential(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> environmentSystems,
    const EmbeddingSettings& settings, std::shared_ptr<PotentialBundle<SCFMode>> activeFockMatrix, bool topDown,
    std::shared_ptr<GridController> supersystemgrid, double gridCutOff,
    std::vector<std::shared_ptr<EnergyComponentController>> allEConts)
  : Potential<SCFMode>(activeSystem->getBasisController()),
    _activeSystem(activeSystem),
    _topDown(topDown),
    _supersystemGrid(supersystemgrid) {
  for (auto e : environmentSystems)
    _environmentSystems.push_back(e);
  activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
      ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);

  _huzinagaProjection = std::make_shared<HuzinagaProjectionPotential<SCFMode>>(
      activeSystem, environmentSystems, settings, activeFockMatrix, topDown, supersystemgrid, gridCutOff, allEConts);

  _s_ABs = SystemSplittingTools<SCFMode>::getProjectedSubsystems(
      activeSystem, environmentSystems, (settings.truncateProjector) ? settings.projecTruncThresh : -10.0);

  _envDensityCont = SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(environmentSystems, topDown);

  EmbeddingSettings envEmbeddingSettings = settings;
  _adjustedSettings_ptr = adjustSettings(envEmbeddingSettings);
  for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
    auto pickedSystem = _environmentSystems[iEnv].lock();
    if (settings.longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE) {
      if (!topDown) {
        LocalizationTask locTask(pickedSystem);
        locTask.run();
      }
      const auto& envCoeff = pickedSystem->template getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          envCoeff, pickedSystem->getOneElectronIntegralController()->getOverlapIntegrals(),
          pickedSystem->getAtomCenteredBasisController()->getBasisIndices());
      auto distantOrbitals = SystemSplittingTools<SCFMode>::selectDistantOrbitals(
          orbitalPopulations, activeSystem, pickedSystem, settings.basisFunctionRatio, settings.borderAtomThreshold);
      auto nonOrthogonalDensityMatrix =
          SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(pickedSystem, distantOrbitals);
      _notProjectedEnvDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(*nonOrthogonalDensityMatrix));
    } // if settings.longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE
  }   // for iEnv
}

template<Options::SCF_MODES SCFMode>
void HoffmannProjectionPotential<SCFMode>::finishSetup() {
  for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
    if (_s_ABs[iEnv]) {
      auto pickedSystem = _environmentSystems[iEnv].lock();
      auto pickedDensityMatrixController = _envDensityCont[iEnv];
      auto otherDensityMatrixController = getOtherDensityMatrixController(
          pickedDensityMatrixController,
          _activeSystem.lock()->getElectronicStructure<SCFMode>()->getDensityMatrixController(), _envDensityCont);
      auto otherSystemController = getOtherSystemController(pickedSystem);
      _embeddingBundle.push_back(FDEPotentialBundleFactory<SCFMode>::produce(
          pickedSystem, pickedDensityMatrixController, otherSystemController, otherDensityMatrixController, _adjustedSettings_ptr,
          (_supersystemGrid) ? _supersystemGrid : pickedSystem->getGridController(), nullptr, _topDown));
    } // if _s_ABs[iEnv]
    else {
      _embeddingBundle.push_back(nullptr);
    }
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<EmbeddingSettings> HoffmannProjectionPotential<SCFMode>::adjustSettings(EmbeddingSettings& settings) {
  settings.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  if (settings.truncateProjector) {
    settings.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  return std::make_shared<EmbeddingSettings>(settings);
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<SystemController>>
HoffmannProjectionPotential<SCFMode>::getOtherSystemController(std::shared_ptr<SystemController> pickedSystem) {
  auto activeSystem = _activeSystem.lock();
  std::vector<std::shared_ptr<SystemController>> otherSystemController = {activeSystem};
  for (auto weak : _environmentSystems) {
    auto env = weak.lock();
    if (env != pickedSystem)
      otherSystemController.push_back(env);
  }
  return otherSystemController;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> HoffmannProjectionPotential<SCFMode>::getOtherDensityMatrixController(
    std::shared_ptr<DensityMatrixController<SCFMode>> pickedDensityMatrixController,
    std::shared_ptr<DensityMatrixController<SCFMode>> activeDensityMatrixController,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> environmentDensityMatrixController) {
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> otherDensityMatrixController = {activeDensityMatrixController};
  for (auto env : environmentDensityMatrixController) {
    if (env != pickedDensityMatrixController)
      otherDensityMatrixController.push_back(env);
  }
  return otherDensityMatrixController;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HoffmannProjectionPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    if (!_finishedSetup) {
      this->finishSetup();
      _finishedSetup = true;
    }
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };

    pot += _huzinagaProjection->getMatrix();

    for (unsigned int iEnv = 0; iEnv < _environmentSystems.size(); ++iEnv) {
      if (_s_ABs[iEnv]) {
        auto environmentSystem = _environmentSystems[iEnv].lock();
        DensityMatrix<SCFMode> projectedMatrix = _envDensityCont[iEnv]->getDensityMatrix();
        if (_notProjectedEnvDensities.size() > 0)
          projectedMatrix -= _notProjectedEnvDensities[iEnv]->getDensityMatrix();
        // A factor of one half for restricted to account for the factor of two in the density matrix.
        double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;
        const FockMatrix<SCFMode> f_BB = _embeddingBundle[iEnv]->getFockMatrix(
            _envDensityCont[iEnv]->getDensityMatrix(), std::make_shared<EnergyComponentController>());
        const Eigen::MatrixXd s_AB = *_s_ABs[iEnv];
        for_spin(f_BB, pot, projectedMatrix) {
          pot_spin +=
              scfFactor * scfFactor * (s_AB * projectedMatrix_spin * f_BB_spin * projectedMatrix_spin * s_AB.transpose());
        };
      } // if _s_ABs[iEnv]
    }   // for iEnv
  }     // if !_potential
  return *_potential;
}
template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HoffmannProjectionPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr += _huzinagaProjection->getGeomGradients();
  return gradientContr;
}

template<Options::SCF_MODES SCFMode>
double HoffmannProjectionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  if (!_potential)
    getMatrix();
  energy += _huzinagaProjection->getEnergy(P);
  return energy;
}

template class HoffmannProjectionPotential<Options::SCF_MODES::RESTRICTED>;
template class HoffmannProjectionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
