/**
 * @file LRSCFAnalysis.h
 *
 * @date Dec. 18, 2018
 * @author Michael Boeckers
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

#ifndef LRSCF_LRSCFANALYSIS
#define LRSCF_LRSCFANALYSIS

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class LRSCFController;

template<Options::SCF_MODES SCFMode>
/**
 * @class LRSCFAnalysis LRSCFAnalysis.h
 * Prints dominant contributions for each determined excitation.
 */
class LRSCFAnalysis {
 public:
  /**
   * @brief Prints dominant contributions
   * @param lrscf The LRSCFController.
   * @param eigenvectors Eigenvectors of the response problem.
   * @param eigenvalues Eigenvalues of the response problem
   * @param contribThresh Threshold to determine which contributions are dominant.
   */
  static void printDominantContributions(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                         const std::vector<Eigen::MatrixXd>& eigenvectors,
                                         const Eigen::VectorXd& eigenvalues, const double contribThresh = 0.9);
};

} /* namespace Serenity */
#endif /* LRSCF_LRSCFANALYSIS */
