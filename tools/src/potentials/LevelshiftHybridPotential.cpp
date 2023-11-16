/**
 * @file   LevelshiftHybridPotential.cpp
 *
 * @date Nov 10, 2018
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
#include "potentials/LevelshiftHybridPotential.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/BasisFunctionProvider.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/SPMatrix.h"
#include "energies/EnergyComponentController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h" //Overlap integrals
#include "misc/SystemSplittingTools.h"
#include "potentials/LevelshiftPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "system/SystemController.h"
#include "tasks/LocalizationTask.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
LevelshiftHybridPotential<SCFMode>::LevelshiftHybridPotential(
    std::shared_ptr<SystemController> activeSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
    const double levelShiftParameter, bool isHybrid, double basisFunctionRatio, double borderAtomThreshold,
    std::shared_ptr<GridController> grid, Functional functional, bool localizedEnv,
    std::pair<bool, std::vector<std::shared_ptr<EnergyComponentController>>> eCon)
  : Potential<SCFMode>(activeSystem->getBasisController()), _isHybrid(isHybrid) {
  /*
   * TODO: If the density of the environment changes for some reason during the lifetime of an instance
   *       of this class, the non-orthogonal orbitals are not selected again. This may lead to problems.
   *       However, this is a rare case which should in general not happen.
   */
  // Set up notifying system.
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  for (const auto& envSys : envSystems) {
    envSys->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
        ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
    DensityMatrix<SCFMode> projectedDensity = envSys->getElectronicStructure<SCFMode>()->getDensityMatrix();
    if (_isHybrid) {
      // Search for distant orbitals
      const auto& envCoeff = envSys->getActiveOrbitalController<SCFMode>()->getCoefficients();
      if (!localizedEnv) {
        LocalizationTask locTask(envSys);
        locTask.run();
      }
      auto orbitalPopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
          envCoeff, envSys->getOneElectronIntegralController()->getOverlapIntegrals(),
          envSys->getAtomCenteredBasisController()->getBasisIndices());
      auto distantOrbitals = SystemSplittingTools<SCFMode>::selectDistantOrbitals(
          orbitalPopulations, activeSystem, envSys, basisFunctionRatio, borderAtomThreshold);
      auto nonOrthogonalDensityMatrix =
          SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(envSys, distantOrbitals);
      projectedDensity -= *nonOrthogonalDensityMatrix;
      _notProjectedEnvDensities.push_back(std::make_shared<DensityMatrixController<SCFMode>>(*nonOrthogonalDensityMatrix));
    } // if isHybrid
    auto projectedDensityMatrixController = std::make_shared<DensityMatrixController<SCFMode>>(projectedDensity);
    // Build level-shift potentials
    _levelshiftPotentials.push_back(std::make_shared<LevelshiftPotential<SCFMode>>(
        activeSystem->getBasisController(), projectedDensityMatrixController, levelShiftParameter));
  } // for envSys
  // Build NAddFuncPotential
  // If false _naddFuncPotential remains nullptr.
  if (_isHybrid) {
    // Fock matrix contribution will now depend on the active density!
    activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
        ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
    _naddFuncPotential = std::make_shared<NAddFuncPotential<SCFMode>>(
        activeSystem, activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        _notProjectedEnvDensities, grid, functional, eCon);
  }
}
template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& LevelshiftHybridPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    for (const auto& levelShift : _levelshiftPotentials)
      *_potential += levelShift->getMatrix();
    if (_naddFuncPotential)
      *_potential += _naddFuncPotential->getMatrix();
  }
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double LevelshiftHybridPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_potential)
    this->getMatrix();
  double energy = 0.0;
  for (const auto& levelShift : _levelshiftPotentials)
    energy += levelShift->getEnergy(P);
  if (_naddFuncPotential)
    energy += _naddFuncPotential->getEnergy(P);
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd LevelshiftHybridPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr.setZero();

  return gradientContr;
}

template class LevelshiftHybridPotential<Options::SCF_MODES::RESTRICTED>;
template class LevelshiftHybridPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
