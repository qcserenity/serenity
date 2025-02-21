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
#include "geometry/Point.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

namespace Options {
enum class LRSCF_TYPE;
} // namespace Options

class Point;

struct LRSCFTaskSettings;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

class SystemController;

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
   * @brief Prints information about the response problem (used functional etc).
   * @param lrscf The LRSCFController.
   * @return The orbital-energy differences, i.e. the leading term of the diagonal of the response matrix.
   */
  static Eigen::VectorXd getDiagonal(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  /**
   * @brief Returns the gauge origin for the given systems and settings.
   * @param settings The LRSCFTaskSettings.
   * @param act The active systems.
   * @param env The environment systems.
   */
  static Point getGaugeOrigin(const LRSCFTaskSettings& settings, const std::vector<std::shared_ptr<SystemController>>& act,
                              const std::vector<std::shared_ptr<SystemController>>& env);

  /**
   * @brief Returns the molecular weight and nuclear dipole moment.
   * @param act The active systems.
   * @param env The environment systems.
   * @param gaugeOrigin The gauge origin.
   * @param molWeight The molecular weight (to be calculated).
   * @param nucDipoleMoment The nuclear dipole moment (to be calculated).
   */
  static void calculateMolecularWeightandNuclearDipoleMoment(const std::vector<std::shared_ptr<SystemController>>& act,
                                                             const std::vector<std::shared_ptr<SystemController>>& env,
                                                             Point gaugeOrigin, double& molWeight,
                                                             Eigen::Vector3d& nucDipoleMoment);

  /**
   * @brief Sets up transformation matrices for the coupled response problem.
   * @param eigenvalues The eigenvalues of the preceeding FDEu tasks are put here.
   * @param settings The LRSCFTaskSettings.
   * @param couplingPatternMatrix The special coupling pattern.
   * @param referenceLoadingType Determine what is being coupled.
   * @param lrscf The LRSCF controller.
   * @param act The active systems.
   * @param nDimension The dimension of the response problem.
   */
  static std::shared_ptr<std::vector<Eigen::MatrixXd>>
  setupFDEcTransformation(Eigen::VectorXd& eigenvalues, const LRSCFTaskSettings& settings,
                          const Eigen::MatrixXi& couplingPatternMatrix,
                          const std::vector<Options::LRSCF_TYPE>& referenceLoadingType,
                          const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                          const std::vector<std::shared_ptr<SystemController>>& act, const unsigned int nDimension);

  /**
   * @brief Sets up all LRSCFController that take part in a response calculation.
   * @param settings The LRSCFTaskSettings.
   * @param act The active systems.
   * @param env The environment systems.
   * @param response The LRSCF controller.
   */
  static void setupLRSCFController(const LRSCFTaskSettings& settings,
                                   const std::vector<std::shared_ptr<SystemController>>& act,
                                   const std::vector<std::shared_ptr<SystemController>>& env,
                                   const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscfAll);

  /**
   * @brief Prepares stability analysis to be done.
   * @param lrscf The LRSCFController.
   * @param settings LRSCFTaskSettings (are modified).
   */
  static void prepareStabilityAnalysis(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                       LRSCFTaskSettings& settings);

  /**
   * @brief Prepares stability analysis to be done.
   * @param lrscf The LRSCFController.
   * @param settings LRSCFTaskSettings (are modified).
   */
  static void printApproximateCoulombInfo(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                          LRSCFTaskSettings& settings);
};

} /* namespace Serenity */

#endif /* LRSCF_LRSCFSETUP */
