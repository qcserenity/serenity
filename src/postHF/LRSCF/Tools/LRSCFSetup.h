/**
 * @file LRSCFSetup.h
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

#ifndef LRSCF_LRSCFSETUP
#define LRSCF_LRSCFSETUP

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"

namespace Serenity {

class Point;

template<Options::SCF_MODES SCFMode>
/**
 * @class LRSCFSetup
 * A static class supposed to keep the task neat and clean.\n
 *
 * For instance, prints information about the response calculation and sets up the\n
 * FDEc transformation matrices upon request.
 */
class LRSCFSetup {
 public:
  /**
   * @brief Prints information about the response problem (used functional etc).
   * @param lrscf The LRSCFController.
   * @param settings The LRSCF settings of the LRSCFTask.
   * @param envSys The systemcontroller of the environment systems.
   * @param type The LRSCF type (Isolated, FDEu, FDEc).
   */
  static void printInfo(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf, const LRSCFTaskSettings& settings,
                        const std::vector<std::shared_ptr<SystemController>>& envSys, const Options::LRSCF_TYPE type);

  /**
   * @brief Returns the gauge origin for the given systems and settings.
   * @param settings The LRSCFTaskSettings.
   * @param act The active systems.
   * @param env The environment systems.
   */
  static Point getGaugeOrigin(const LRSCFTaskSettings& settings, const std::vector<std::shared_ptr<SystemController>>& act,
                              const std::vector<std::shared_ptr<SystemController>>& env);

  /**
   * @brief Sets up transformation matrices for the coupled response problem.
   * @param couplingPatternMatrix The special coupling pattern.
   * @param referenceLoadingType Determine what is being coupled.
   * @param lrscf The LRSCF controller.
   * @param act The active systems.
   * @param nDimension The dimension of the response problem.
   */
  static shared_ptr<std::vector<Eigen::MatrixXd>>
  setupFDEcTransformation(const LRSCFTaskSettings& settings, const Eigen::MatrixXi& couplingPatternMatrix,
                          const std::vector<Options::LRSCF_TYPE>& referenceLoadingType,
                          const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                          const std::vector<std::shared_ptr<SystemController>>& act, const unsigned int nDimension);

  /**
   * @brief Sets up all LRSCFController that take part in a response calculation.
   * @param settings The LRSCFTaskSettings.
   * @param couplingPatternMatrix The special coupling pattern.
   * @param act The active systems.
   * @param env The environment systems.
   * @param response The LRSCF controller.
   * @param type The LRSCF type (Isolated, FDEu, FDEc).
   */
  static void setupLRSCFController(const LRSCFTaskSettings& settings, const Eigen::MatrixXi& couplingPatternMatrix,
                                   const std::vector<std::shared_ptr<SystemController>>& act,
                                   const std::vector<std::shared_ptr<SystemController>>& env,
                                   const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& response,
                                   const Options::LRSCF_TYPE type);
};

} /* namespace Serenity */

#endif /* LRSCF_LRSCFSETUP */
