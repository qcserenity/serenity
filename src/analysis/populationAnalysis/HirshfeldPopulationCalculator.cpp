/**
 * @file HirshfeldPopulationCalculator.cpp
 *
 * @date Mar 11, 2016
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
#include "analysis/populationAnalysis/HirshfeldPopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "geometry/Geometry.h"
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HirshfeldPopulationCalculator<SCFMode>::HirshfeldPopulationCalculator(std::shared_ptr<SystemController> systemController,
                                                                      std::shared_ptr<DensityOnGridController<SCFMode>> densOnGrid)
  : _system(systemController), _densOnGrid(densOnGrid), _atomPopulations(nullptr) {
  this->_densOnGrid->addSensitiveObject(ObjectSensitiveClass<DensityOnGrid<SCFMode>>::_self);
}

template<Options::SCF_MODES SCFMode>
void HirshfeldPopulationCalculator<SCFMode>::calculateHirshfeldAtomPopulations() {
  // Get grid representation of relaxed molecular density
  auto basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _system->getSettings(), _system->getAtomCenteredBasisController(), _system->getGridController());
  auto densOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
      basisFunctionOnGridController, _system->getSettings().grid.blockAveThreshold);
  GridData<RESTRICTED> convergedSystem = _densOnGrid->getDensityOnGrid().total();

  // Grid weights, points
  auto& weightOfGrid = _system->getGridController()->getWeights();
  auto nBlocks = basisFunctionOnGridController->getNBlocks();
  auto nGridPoints = basisFunctionOnGridController->getNGridPoints();

  auto targetBasisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _system->getSettings(), _system->getBasisController(), _system->getGridController());
  ScalarOperatorToMatrixAdder<RESTRICTED> gridToMat(targetBasisFunctionOnGridController,
                                                    _system->getSettings().grid.blockAveThreshold);

  // Guess using atom densities which are calculated in place using the settings of the target system
  AtomicDensityGuessCalculator guess(GUESSMODES::SCF_INPLACE);
  DensityMatrix<RESTRICTED> proMolecularDensityMat(_system->getAtomCenteredBasisController());
  proMolecularDensityMat += (Eigen::MatrixXd)*guess.calculateInitialDensity(_system, false, false);

  // Indices of basis functions centered on corresponding atoms
  const auto& indices = _system->getAtomCenteredBasisController()->getBasisIndices();

  // Atoms
  auto atoms = _system->getGeometry()->getAtoms();
  auto nAtoms = _system->getGeometry()->getNAtoms();

  // Initialize
  DensityOnGrid<RESTRICTED> proMolDens(_system->getGridController());
  std::vector<DensityOnGrid<RESTRICTED>> atomDensities;
  _atomPopulations.reset(new Eigen::VectorXd(Eigen::VectorXd::Zero(nAtoms)));
  auto& atomPopulations = *_atomPopulations;

  // Create promolecular density and free atom densities
  for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
    std::vector<std::shared_ptr<Atom>> atomVec = {atoms[iAtom]};
    auto atomGeom = std::make_shared<Geometry>(atomVec);
    auto atomBasisController = AtomCenteredBasisControllerFactory::produce(
        atomGeom, _system->getSettings().basis.basisLibPath, _system->getSettings().basis.makeSphericalBasis, false,
        _system->getSettings().basis.firstECP, "");
    auto atomBasisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        _system->getSettings(), atomBasisController, _system->getGridController());
    auto atomDensOnGridCalculator = std::make_shared<DensityOnGridCalculator<RESTRICTED>>(
        atomBasisFunctionOnGridController, _system->getSettings().grid.blockAveThreshold);

    /* Extracting atom density from promolecular density*/
    DensityMatrix<RESTRICTED> atomDensityMatrix(atomBasisController);
    atomDensityMatrix = proMolecularDensityMat.block(indices[iAtom].first, indices[iAtom].first,
                                                     indices[iAtom].second - indices[iAtom].first,
                                                     indices[iAtom].second - indices[iAtom].first);

    // Get grid representation
    atomDensities.push_back(atomDensOnGridCalculator->calcDensityOnGrid(atomDensityMatrix));
    proMolDens += atomDensities[iAtom];
  }

  for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
    // Prepare variables
    Eigen::VectorXd thisAtomPopulation = Eigen::VectorXd::Zero(omp_get_max_threads());
    // Calculate
#pragma omp for schedule(static)
    for (unsigned int block = 0; block < nBlocks; block++) {
      unsigned int blockEnd;
      if (block == (nBlocks - 1)) {
        blockEnd = nGridPoints;
      }
      else {
        blockEnd = basisFunctionOnGridController->getFirstIndexOfBlock(block + 1);
      }
      auto blockStart = basisFunctionOnGridController->getFirstIndexOfBlock(block);
      for (unsigned int i = blockStart; i < blockEnd; ++i) {
        if (proMolDens[i] < 1.0e-9)
          continue;
        const double weightPerAtomGrid = atomDensities[iAtom][i] / proMolDens[i];
        thisAtomPopulation[omp_get_thread_num()] += convergedSystem[i] * weightPerAtomGrid * weightOfGrid[i];
      }
    }

    atomPopulations[iAtom] = thisAtomPopulation.sum();

  } // iAtom
}
template class HirshfeldPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class HirshfeldPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity
