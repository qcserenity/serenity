/**
 * @file BeckePopulationCalculator.cpp
 *
 * @date Jun 16, 2020
 * @author: Patrick Eschenbach
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "analysis/populationAnalysis/BeckePopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGrid.h"
#include "grid/AtomCenteredGridController.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "grid/GridController.h"
#include "settings/Options.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
BeckePopulationCalculator<SCFMode>::BeckePopulationCalculator(std::shared_ptr<SystemController> systemController,
                                                              std::shared_ptr<DensityMatrix<SCFMode>> density)
  : _system(systemController), _density(density) {
  assert(_system);
  assert(_density);
}

template<Options::SCF_MODES SCFMode>
void BeckePopulationCalculator<SCFMode>::calculateBeckeAtomPopulations() {
  Settings beckeSettings = _system->getSettings();
  beckeSettings.grid.gridType = Options::GRID_TYPES::BECKE;
  beckeSettings.grid.radialGridType = Options::RADIAL_GRID_TYPES::BECKE;
  beckeSettings.grid.accuracy = 7;
  beckeSettings.grid.smallGridAccuracy = 7;
  beckeSettings.grid.gridPointSorting = false;
  auto atomGridController = AtomCenteredGridControllerFactory::produce(_system->getGeometry(), beckeSettings.grid);
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(beckeSettings.grid.blocksize, beckeSettings.grid.basFuncRadialThreshold,
                                                    0, _system->getBasisController(), atomGridController);
  auto densOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(basFuncOnGridController, 0.0);
  auto points = atomGridController->getGridPoints();
  auto atomWeights = atomGridController->getWeights();
  auto atomIndices = atomGridController->getAtomGrid()->getGridIndicesOfAtoms();
  unsigned nAtoms = _system->getGeometry()->getNAtoms();
  auto densMatController = std::make_shared<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>((*_density).total());
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>(
      densOnGridCalculator, densMatController);
  auto densOnGrid = DensityOnGrid<Options::SCF_MODES::RESTRICTED>(atomGridController);
  densOnGrid = densOnGridController->getDensityOnGrid();
  auto atomPopulations = std::make_shared<Eigen::VectorXd>(Eigen::VectorXd::Zero(nAtoms));
#pragma omp parallel for schedule(dynamic)
  for (unsigned iAtom = 0; iAtom < nAtoms; ++iAtom) {
    double thisAtomPopulation = 0.0;
    for (unsigned iPoint = atomIndices[iAtom].first; iPoint < atomIndices[iAtom].second; ++iPoint) {
      thisAtomPopulation += densOnGrid[iPoint] * atomWeights[iPoint];
    }
    (*atomPopulations)[iAtom] = thisAtomPopulation;
  }
  _atomPopulations = atomPopulations;
}

template<Options::SCF_MODES SCFMode>
void BeckePopulationCalculator<SCFMode>::calculateBeckeSpinPopulations() {
  Settings beckeSettings = _system->getSettings();
  beckeSettings.grid.gridType = Options::GRID_TYPES::BECKE;
  beckeSettings.grid.radialGridType = Options::RADIAL_GRID_TYPES::BECKE;
  beckeSettings.grid.accuracy = 7;
  beckeSettings.grid.smallGridAccuracy = 7;
  beckeSettings.grid.gridPointSorting = false;
  auto atomGridController = AtomCenteredGridControllerFactory::produce(_system->getGeometry(), beckeSettings.grid);
  auto basFuncOnGridController =
      BasisFunctionOnGridControllerFactory::produce(beckeSettings.grid.blocksize, beckeSettings.grid.basFuncRadialThreshold,
                                                    0, _system->getBasisController(), atomGridController);
  auto densOnGridCalculator =
      std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(basFuncOnGridController, 0.0);
  auto points = atomGridController->getGridPoints();
  auto atomWeights = atomGridController->getWeights();
  auto atomIndices = atomGridController->getAtomGrid()->getGridIndicesOfAtoms();
  unsigned nAtoms = _system->getGeometry()->getNAtoms();
  auto densMatController =
      std::make_shared<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>((*_density).difference());
  auto densOnGridController = std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>(
      densOnGridCalculator, densMatController);
  auto densOnGrid = DensityOnGrid<Options::SCF_MODES::RESTRICTED>(atomGridController);
  densOnGrid = densOnGridController->getDensityOnGrid();
  auto spinPopulations = std::make_shared<Eigen::VectorXd>(Eigen::VectorXd::Zero(nAtoms));
#pragma omp parallel for schedule(dynamic)
  for (unsigned iAtom = 0; iAtom < nAtoms; ++iAtom) {
    double thisAtomPopulation = 0.0;
    for (unsigned iPoint = atomIndices[iAtom].first; iPoint < atomIndices[iAtom].second; ++iPoint) {
      thisAtomPopulation += densOnGrid[iPoint] * atomWeights[iPoint];
    }
    (*spinPopulations)[iAtom] = thisAtomPopulation;
  }
  _spinPopulations = spinPopulations;
}
template class BeckePopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class BeckePopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity