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
using namespace std;

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
    const CoefficientMatrix<SCFMode>& coefficients, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
    const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices) {
  assert(isDefinedInSameBasis(coefficients, overlapMatrix));
  /*
   * Get in data locally
   */
  const unsigned int nAtoms = atomBasisIndices.size();
  const unsigned int nBasisFunctions = overlapMatrix.getNBasisFunctions();
  /* Prepare output */
  SPMatrix<SCFMode> atomwiseOrbitalPopulations(nAtoms, nBasisFunctions);

  // Does not HAVE to be equal to nBasisFunctions in principle but is assumed in many places.
  unsigned int nOrbitals;
  for_spin(coefficients) {
    nOrbitals = coefficients_spin.cols();
  };
  for (unsigned int i = 0; i < nOrbitals; ++i) {
    for_spin(atomwiseOrbitalPopulations, coefficients) {
      const auto basisFunctionOrbitalPopulations = calculateOrbitalPopulations(coefficients_spin.col(i), overlapMatrix);
      /*
       * Loop over atoms
       */
      for (unsigned int k = 0; k < nAtoms; ++k) {
        // Find out which matrix entries belong to basis functions of this atom
        atomwiseOrbitalPopulations_spin(k, i) = 0.0;
        for (unsigned int mu = atomBasisIndices[k].first; mu < atomBasisIndices[k].second; ++mu) {
          atomwiseOrbitalPopulations_spin(k, i) += basisFunctionOrbitalPopulations[mu];
        }
      }
    };
  }
  return atomwiseOrbitalPopulations;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd MullikenPopulationCalculator<SCFMode>::calculateOrbitalPopulations(
    const Eigen::VectorXd& orbitalcoeffitients, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix) {
  assert((int)orbitalcoeffitients.size() == overlapMatrix.rows());
  /*
   * Get in data locally
   */
  const unsigned int nBasisFunctions = overlapMatrix.getNBasisFunctions();
  /* Prepare output */
  Eigen::VectorXd populations = Eigen::VectorXd::Zero(nBasisFunctions);

  /*
   * Loop over basis functions
   */
  for (unsigned int i = 0; i < nBasisFunctions; ++i) {
    for (unsigned int j = 0; j < nBasisFunctions; ++j) {
      populations[i] += orbitalcoeffitients(i) * orbitalcoeffitients(j) * overlapMatrix(i, j);
    }
  }
  return populations;
}

template class MullikenPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class MullikenPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
