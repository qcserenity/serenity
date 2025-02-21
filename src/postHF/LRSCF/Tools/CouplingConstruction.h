/**
 * @file CouplingConstruction.h
 *
 * @date Jan. 10, 2019
 * @author Johannes Tölle
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

#ifndef LRSCF_COUPLINGCONSTRUCTION
#define LRSCF_COUPLINGCONSTRUCTION

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Tools/SigmaCalculator.h"
#include "settings/ElectronicStructureOptions.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

namespace Options {
enum class LRSCF_TYPE;
} // namespace Options

struct LRSCFTaskSettings;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class CouplingConstruction
 * A static class supposed to keep the task neat and clean.\n
 *
 * Performs the partial response matrix construction to save\n
 * some time during a one-step FDEc procedure.\n
 * This routine cannot be applied for full FDEc!
 */
template<Options::SCF_MODES SCFMode>
class CouplingConstruction {
 public:
  /**
   * @brief Performs the partial response-matrix construction.
   * @param lrscf The LRSCFController.
   * @param settings The LRSCF settings of the LRSCFTask.
   * @param sigmaCalculator The sigmacalculator to form the partial matrix-vector product.
   * @param eigenvectors The eigenvectors to hold the solution.
   * @param eigenvalues The eigenvalues to hold the solution.
   */
  static void solve(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, const LRSCFTaskSettings& settings,
                    std::vector<Options::LRSCF_TYPE> referenceLoadingType, SigmaCalculator sigmaCalculator,
                    std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors, Eigen::VectorXd& eigenvalues);
};

} /* namespace Serenity */

#endif /* LRSCF_COUPLINGCONSTRUCTION */
