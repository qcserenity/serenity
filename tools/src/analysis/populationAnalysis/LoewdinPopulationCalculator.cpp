/**
 * @file   LoewdinPopulationCalculator.cpp
 *
 * @date   Oct 05, 2021
 * @author Niklas Niemeyer
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
#include "analysis/populationAnalysis/LoewdinPopulationCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "system/SystemController.h"
/* Include Std and External Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
LoewdinPopulationCalculator<SCFMode>::calculateLoewdinPopulations(std::shared_ptr<SystemController> systemController) {
  return calculateAtomPopulations(systemController->getElectronicStructure<SCFMode>()->getDensityMatrix(),
                                  systemController->getOneElectronIntegralController()->getOverlapIntegrals(),
                                  systemController->getAtomCenteredBasisController()->getBasisIndices());
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> LoewdinPopulationCalculator<SCFMode>::calculateAtomPopulations(
    const DensityMatrix<SCFMode>& densityMatrix, const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
    const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices) {
  assert(isDefinedInSameBasis(densityMatrix, overlapMatrix));

  Eigen::MatrixXd sqrtS = mSqrt_Sym(overlapMatrix);

  const unsigned int nAtoms = atomBasisIndices.size();
  const auto basisFunctionPopulations = calculateBasisFunctionPopulations(densityMatrix, sqrtS);
  SpinPolarizedData<SCFMode, Eigen::VectorXd> atomPopulations(Eigen::VectorXd::Zero(nAtoms));

  for_spin(atomPopulations, basisFunctionPopulations) {
    for (unsigned int i = 0; i < nAtoms; ++i) {
      for (unsigned int j = atomBasisIndices[i].first; j < atomBasisIndices[i].second; ++j) {
        atomPopulations_spin[i] += basisFunctionPopulations_spin[j];
      }
    }
  };
  return atomPopulations;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
LoewdinPopulationCalculator<SCFMode>::calculateBasisFunctionPopulations(const DensityMatrix<SCFMode>& densityMatrix,
                                                                        const Eigen::MatrixXd& sqrtS) {
  const unsigned int nBasisFunctions = densityMatrix.getNBasisFunctions();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> populations(Eigen::VectorXd::Zero(nBasisFunctions));

  for_spin(populations, densityMatrix) {
    Eigen::MatrixXd tmp = sqrtS * densityMatrix_spin * sqrtS;
    for (unsigned int i = 0; i < nBasisFunctions; ++i) {
      populations_spin[i] = tmp(i, i);
    }
  };
  return populations;
}

template class LoewdinPopulationCalculator<Options::SCF_MODES::RESTRICTED>;
template class LoewdinPopulationCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
