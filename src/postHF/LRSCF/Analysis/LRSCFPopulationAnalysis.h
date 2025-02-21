/**
 * @file LRSCFPopulationAnalysis.h
 * @author Niklas Niemeyer, Anton Rikus
 *
 * @date October 10, 2021
 *
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

#ifndef LRSCFPOPULATIONANALYSIS_H_
#define LRSCFPOPULATIONANALYSIS_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class LRSCFController;
/**
 * @class LRSCFPopulationAnalysis
 * @brief Calculates transition, particle and hole populations and stores them on disk.
 */
template<Options::SCF_MODES SCFMode>
class LRSCFPopulationAnalysis {
 public:
  /**
   * @brief Transforms transition density matrices to the AO basis, performs a Löwdin population analysis
   * of it and stores it on disk for later analysis.
   *
   * Also calculates hole and particle density matrices and transforms them to the AO basis to then perform
   * a Löwdin population analysis of it. Outputs particle and hole charges and stores the correlation matrix on disk.
   * The correlation matrix's order on disk is in such a way that the probability that a particle state is at atom 1
   * while the hole state is at atom 1 as well is in the lower left corner. Particle states are on the y-axis while hole
   * states are on the x-axis.
   *
   * @param lrscf The vector of LRSCF controllers storing all information.
   * @param densityMatrices The vector of transition density matrices (simply the excitation vectors in TDDFT).
   */
  static void calculateTransitionCharges(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                         std::vector<Eigen::MatrixXd>& densityMatrices);
};

} /* namespace Serenity */
#endif /* LRSCFPOPULATIONANALYSIS_H_ */
