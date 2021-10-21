/**
 * @file   MullikenPopulationCalculator.cpp
 *
 * @date   Mar 11, 2014
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
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "data/matrices/SPMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/Matrix.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <cassert>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
MullikenPopulationCalculator<SCFMode>::calculateMullikenPopulations(std::shared_ptr<SystemController> systemController) {
  return calculateAtomPopulations(systemController->getElectronicStructure<SCFMode>()->getDensityMatrix(),
                                  systemController->getOneElectronIntegralController()->getOverlapIntegrals(),
                                  systemController->getAtomCenteredBasisController()->getBasisIndices());
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> MullikenPopulationCalculator<SCFMode>::calculateAtomPopulations(
    const DensityMatrix<SCFMode>& densityMatrix, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
    const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices) {
  assert(isDefinedInSameBasis(densityMatrix, overlapMatrix));
  /*
   * Get in data locally
   */
  const unsigned int nAtoms = atomBasisIndices.size();
  const auto basisFunctionPopulations = calculateBasisFunctionPopulations(densityMatrix, overlapMatrix);
  /* Prepare output */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> atomPopulations(Eigen::VectorXd::Zero(nAtoms));
  /*
   * Loop over atoms
   */
  for_spin(atomPopulations, basisFunctionPopulations) {
    for (unsigned int i = 0; i < nAtoms; ++i) {
      // Find out which matrix entries belong to basis functions of this atom
      for (unsigned int j = atomBasisIndices[i].first; j < atomBasisIndices[i].second; ++j) {
        atomPopulations_spin[i] += basisFunctionPopulations_spin[j];
      }
    }
  };
  return atomPopulations;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> MullikenPopulationCalculator<SCFMode>::calculateBasisFunctionPopulations(
    const DensityMatrix<SCFMode>& densityMatrix, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix) {
  assert(isDefinedInSameBasis(densityMatrix, overlapMatrix));
  /*
   * Get in data locally
   */
  const unsigned int nBasisFunctions = densityMatrix.getNBasisFunctions();
  /* Prepare output */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> populations(Eigen::VectorXd::Zero(nBasisFunctions));
  /*
   * Loop over basis functions
   */
  for_spin(populations, densityMatrix) {
    for (unsigned int i = 0; i < nBasisFunctions; ++i) {
      for (unsigned int j = 0; j < nBasisFunctions; ++j) {
        populations_spin[i] += densityMatrix_spin(i, j) * overlapMatrix(i, j);
      }
    }
  };
  return populations;
}
template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
    const SPMatrix<SCFMode>& coefficients, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
    const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices) {
  // Resize later!
  SPMatrix<SCFMode> atomwiseOrbitalPopulations(1, 1);
  const unsigned int nAtoms = atomBasisIndices.size();
  for_spin(coefficients, atomwiseOrbitalPopulations) {
    const unsigned int nOrbs = coefficients.cols();
    atomwiseOrbitalPopulations_spin.resize(nAtoms, coefficients.cols());
    // For each MO: Calculate basis function population and sum them up according to the atoms.
    Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iOrb = 0; iOrb < nOrbs; ++iOrb) {
      const Eigen::VectorXd basisFunctionOrbitalPopulations =
          calculateOrbitalPopulations(coefficients_spin.col(iOrb), overlapMatrix);
      for (unsigned int k = 0; k < nAtoms; ++k) {
        // Find out which matrix entries belong to basis functions of this atom and sum them up.
        unsigned nBasOnAtom = atomBasisIndices[k].second - atomBasisIndices[k].first;
        atomwiseOrbitalPopulations_spin(k, iOrb) =
            basisFunctionOrbitalPopulations.segment(atomBasisIndices[k].first, nBasOnAtom).sum();
      }
    }
    Eigen::setNbThreads(0);
  };
  return atomwiseOrbitalPopulations;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd MullikenPopulationCalculator<SCFMode>::calculateOrbitalPopulations(
    const Eigen::VectorXd& orbitalcoeffitients, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix) {
  Eigen::VectorXd populations =
      ((orbitalcoeffitients * orbitalcoeffitients.transpose()).array() * overlapMatrix.array()).rowwise().sum();
  return populations;
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>
MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(std::shared_ptr<SystemController> system) {
  return calculateAtomwiseOrbitalPopulations(system->template getActiveOrbitalController<SCFMode>()->getCoefficients(),
                                             system->getOneElectronIntegralController()->getOverlapIntegrals(),
                                             system->getAtomCenteredBasisController()->getBasisIndices());
}

template class MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
