/**
 * @file   GeometryOptimizationTask.h
 *
 * @date   May 21, 2015
 * @author Jan Unsleber
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
#ifndef CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_
#define CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_
/* Include Serenity Internal Headers */
#include "misc/Timing.h"
#include "settings/DFTOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/GeometryOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Optimizer;
/**
 * @class SQNMSettings GeometryOptimizationTask.h
 * @brief Settings for the SQNM algorithm.\n
 *
 *
 * Settings:\n
 * -historyLength    - Maximum number of coordinate and gradient information of
 *                     previous cycles used in the optimization.\n
 * -epsilon          - Threshold to determine significant eigenvalues of the displacement overlap.\n
 * -alpha            - Initial step length to be updated during the optimization.\n
 * -energyThreshold  - Energy threshold to determine whether an optimization step is accepted.\n
 * -trustRadius      - Maximum step length - step will be  scaled to this length if the trust radius is exceeded.\n
 */
using namespace Serenity::Reflection;
struct SQNMSettings {
  SQNMSettings() : historyLength(10), epsilon(1.0e-4), alpha(1.0), energyThreshold(1.0e-6), trustRadius(0.1) {
  }

 public:
  REFLECTABLE((unsigned int)historyLength, (double)epsilon, (double)alpha, (double)energyThreshold, (double)trustRadius)

  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("SQNM")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};

using namespace Serenity::Reflection;
struct GeometryOptimizationTaskSettings {
  GeometryOptimizationTaskSettings()
    : gradType(Options::GRADIENT_TYPES::ANALYTICAL),
      maxCycles(100),
      rmsgradThresh(1.0e-4),
      energyChangeThresh(5.0e-6),
      maxGradThresh(3.0e-4),
      stepThresh(2.0e-3),
      maxStepThresh(4.0e-3),
      numGradStepSize(0.001),
      transInvar(false),
      FaTmaxCycles(50),
      FaTConvThresh(1.0e-6),
      FaTgridCutOff(-1.0),
      optAlgorithm(Options::OPTIMIZATION_ALGORITHMS::BFGS) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  REFLECTABLE((Options::GRADIENT_TYPES)gradType, (unsigned int)maxCycles, (double)rmsgradThresh,
              (double)energyChangeThresh, (double)maxGradThresh, (double)stepThresh, (double)maxStepThresh,
              (double)numGradStepSize, (bool)transInvar, (unsigned int)FaTmaxCycles, (double)FaTConvThresh,
              (double)FaTgridCutOff, (Options::OPTIMIZATION_ALGORITHMS)optAlgorithm)
 public:
  EmbeddingSettings embedding;
  SQNMSettings sqnm;
};
/**
 * @class  GeometryOptimizationTask GeometryOptimizationTask.h
 * @brief  Task to optimize the geometry of the system
 */
template<Options::SCF_MODES SCFMode>
class GeometryOptimizationTask : public Task {
 public:
  /**
   * @param systemController
   */

  GeometryOptimizationTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                           const std::vector<std::shared_ptr<SystemController>>& passiveSystems =
                               std::vector<std::shared_ptr<SystemController>>());
  /**
   * @brief Default destructor.
   */
  virtual ~GeometryOptimizationTask() = default;
  /**
   * @see Task
   */
  void run() override;
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(GeometryOptimizationTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname)) {
      return;
    }
    if (c.sqnm.visitAsBlockSettings(v, blockname)) {
      return;
    }
    throw SerenityError((std::string) "Unknown block in GeometryOptimizationTaskSettings: " + blockname);
  }
  /**
   * @brief The settings/keywords for the GeometryOptimizationTask:
   *        - gradType: The type of the gradient. The following options are available:
   *          - NUMERICAL
   *          - ANALYTICAL (default)
   *        - maxCycles: Maximum number of steps (default: 100)
   *        - gradNormThresh: The gradient convergence threshold (default: 1.0e-5)
   *        - numGradStepSize: The displacement step size used for construction of the numerical gradients (default 0.05)
   *        - transInvar Make gradients translationally invariant.
   *        - FaTmaxCycles The maximum number of FaT iterations.
   *        - FaTConvThresh Convergence threshold for the FaT calculation.
   *        - FaTgridCutOff A distance cutoff for the integration grid used to calculate
   *                        the non-additive energy functional potentials.
   *        -optAlgorithm: The type of optimization algorithm to be used. The following options are available:
   *          - BFGS (default)
   *          - SQNM
   *
   */
  GeometryOptimizationTaskSettings settings;

 private:
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  std::vector<std::shared_ptr<SystemController>> _passiveSystems;
  bool _minimumPrint = false;
  timespec _time;
};

} /* namespace Serenity */

#endif /* CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_ */
