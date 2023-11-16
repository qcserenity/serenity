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
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/DensityOnGridController.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/grid/SupersystemDensityOnGridController.h"
#include "data/matrices/DensityMatrixController.h"
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
  const GridData<SCFMode>& convergedSystem = _densOnGrid->getDensityOnGrid();

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
  proMolecularDensityMat += (Eigen::MatrixXd)*guess.calculateInitialDensity(_system);

  // Indices of basis functions centered on corresponding atoms
  const auto& indices = _system->getAtomCenteredBasisController()->getBasisIndices();

  // Atoms
  auto atoms = _system->getGeometry()->getAtoms();
  auto nAtoms = _system->getGeometry()->getNAtoms();

  // Initialize
  std::vector<std::shared_ptr<DensityOnGridController<RESTRICTED>>> atomDensities;
  std::vector<Eigen::VectorXi> nonNegligibleBlockAtomWise;
  _atomPopulations = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(Eigen::VectorXd::Zero(nAtoms));
  SpinPolarizedData<SCFMode, Eigen::VectorXd>& atomPopulations = *_atomPopulations;

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
    // MB: I will use the density matrix controller here in order to retain the prescreening information
    //     from the BasisFunctionOnGridController.
    auto atomDensityMatrixController = std::make_shared<DensityMatrixController<RESTRICTED>>(atomDensityMatrix);
    auto atomDensityMatrixDensityOnGridController = std::make_shared<DensityMatrixDensityOnGridController<RESTRICTED>>(
        atomDensOnGridCalculator, atomDensityMatrixController, 0);
    atomDensities.push_back(atomDensityMatrixDensityOnGridController);
    atomDensityMatrixDensityOnGridController->getDensityOnGrid();
    nonNegligibleBlockAtomWise.push_back(atomDensityMatrixDensityOnGridController->getNonNegligibleBlocks());
  }
  // Build the pro molecular density from the superposition of the atom densities. We will recycle the FDE-framework
  // here.
  auto proMolDensController = std::make_shared<SupersystemDensityOnGridController<RESTRICTED>>(atomDensities);
  const DensityOnGrid<RESTRICTED>& proMolDens = proMolDensController->getDensityOnGrid();

  for_spin(convergedSystem, atomPopulations) {
    Eigen::MatrixXd threadWiseAtomPops = Eigen::MatrixXd::Zero(nAtoms, omp_get_num_threads());
    // The outer loop is over the grid points block in order to allow the precalculation of
    // the product of density and grid weights as well as efficient vector operations.
#pragma omp for schedule(dynamic)
    for (unsigned int block = 0; block < nBlocks; block++) {
      unsigned int blockEnd;
      if (block == (nBlocks - 1)) {
        blockEnd = nGridPoints;
      }
      else {
        blockEnd = basisFunctionOnGridController->getFirstIndexOfBlock(block + 1);
      }
      auto blockStart = basisFunctionOnGridController->getFirstIndexOfBlock(block);
      const unsigned int blockSize = blockEnd - blockStart;

      // Calculate product of density and grid weights. These will be used for every atom.
      const Eigen::VectorXd weightedDensity = convergedSystem_spin.segment(blockStart, blockSize).array() *
                                              weightOfGrid.segment(blockStart, blockSize).array();
      // Calculate the reciprocal of the pro molecular density.
      // If the pro molecular density is small, we will simply assign its reciprocal zero.
      // Note that the pro molecular density is always equal or larger than the atom densities. Thus, this
      // does not introduce any error.
      Eigen::VectorXd reciProMolDens = Eigen::VectorXd::Zero(blockSize);
      for (unsigned int i = 0; i < blockSize; ++i) {
        const double proMolDens_i = proMolDens(i + blockStart);
        if (proMolDens_i > 1.0e-9)
          reciProMolDens(i) = 1.0 / proMolDens_i;
      }
      // Calculate atom-specific weights and contract with weighted density.
      for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
        // Skip the atom, if there are no non-zero entries in the atom density in this block.
        if (not nonNegligibleBlockAtomWise[iAtom](block))
          continue;
        const Eigen::VectorXd weightPerAtomGrid =
            atomDensities[iAtom]->getDensityOnGrid().segment(blockStart, blockSize).array() * reciProMolDens.array();
        threadWiseAtomPops(iAtom, omp_get_thread_num()) += (weightPerAtomGrid.array() * weightedDensity.array()).sum();
      } // iAtom
    }
    atomPopulations_spin = threadWiseAtomPops.rowwise().sum();
  };
}
template class HirshfeldPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class HirshfeldPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity
