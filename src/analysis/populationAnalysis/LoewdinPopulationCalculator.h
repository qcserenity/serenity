/**
 * @file   LoewdinPopulationCalculator.h
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
#ifndef LOEWDINPOPULATIONCALCULATOR_H_
#define LOEWDINPOPULATIONCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */

class SystemController;
/**
 * @class  LoewdinPopulationCalculator LoewdinPopulationCalculator.h
 * @brief  Performs a population analysis according to Loewdin.
 *
 * I.e. electrons are counted to an atom if they reside in basis functions which are centered
 * on that atom.
 */
template<Options::SCF_MODES SCFMode>
class LoewdinPopulationCalculator {
 public:
  /**
   * @brief Calculates atom-wise atom populations.
   * @param   systemController The system for which atom-wise Löwdin populations are to be calculated.
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd> calculateLoewdinPopulations(std::shared_ptr<SystemController> systemController);

  /**
   * @brief Calculates atom-wise atom populations.
   *
   * @param densityMatrix   The AO density matrix.
   * @param overlapMatrix   The AO overlap matrix (for which the sqrt is computed).
   * @param atomBasisIndices Used to map from basis functions to atoms.
   *                         For each atom there is a first index and an end index
   * @returns The Loewdin population on the atoms defined by the atomBasisIndices.
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd>
  calculateAtomPopulations(const DensityMatrix<SCFMode>& densityMatrix,
                           const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
                           const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices);

  /**
   * @brief Calculates the population of each basis function.
   *
   * @param densityMatrix   The AO density matrix.
   * @param overlapMatrix   The AO overlap matrix (for which the sqrt is computed).
   * @returns The population of each basis function, i.e.
   *          \f$ {\rm Vector}_{\mu} = \sum_{\nu} P_{\mu,\nu} \cdot S_{\mu,\nu} \f$
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd>
  calculateBasisFunctionPopulations(const DensityMatrix<SCFMode>& densityMatrix, const Eigen::MatrixXd& overlapMatrix);
};

} /* namespace Serenity */

#endif /* LOEWDINPOPULATIONCALCULATOR_H_ */
