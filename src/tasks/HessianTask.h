/**
 * @file HessianTask.h
 *
 * @date Feb 7, 2017
 * @author Kevin Klahr
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

#ifndef TASKS_HESSIANTASK_H_
#define TASKS_HESSIANTASK_H_

/* Include Serenity Internal Headers */
#include "settings/DFTOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/GeometryOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;
struct HessianTaskSettings {
  HessianTaskSettings()
    : gradType(Options::GRADIENT_TYPES::ANALYTICAL),
      hessType(Options::HESSIAN_TYPES::NUMERICAL),
      numHessStepSize(0.001),
      numGradStepSize(0.001),
      printToFile(true),
      FaTmaxCycles(50),
      FaTenergyConvThresh(1.0e-6),
      FaTgridCutOff(-1.0) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  REFLECTABLE((Options::GRADIENT_TYPES)gradType, (Options::HESSIAN_TYPES)hessType, (double)numHessStepSize,
              (double)numGradStepSize, (bool)printToFile, (int)FaTmaxCycles, (double)FaTenergyConvThresh, (double)FaTgridCutOff)
 public:
  EmbeddingSettings embedding;
};
/**
 * @class  HessianTask HessianTask.h
 * @brief  Task to calculate cartesian Hessian of the system
 */
template<Options::SCF_MODES SCFMode>
class HessianTask : public Task {
 public:
  /**
   * @param systemController         The system for which the Hessian shall be calculated
   */
  HessianTask(std::vector<std::shared_ptr<SystemController>> activeSystems,
              std::vector<std::shared_ptr<SystemController>> passiveSystems = std::vector<std::shared_ptr<SystemController>>());
  /**
   * @brief Default destructor.
   */
  virtual ~HessianTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(HessianTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in FreezeAndThawTaskSettings: " + blockname);
  }
  /**
   * @brief The settings/keywords for HessianTask:
   *        - gradType: The type of the gradient. The following options are available:
   *          - NUMERICAL
   *          - ANALYTICAL (default)
   *        - hessType: The type of the hessian. The following options are available:
   *          - NUMERICAL (default)
   *          - ANALYTICAL
   *        - numGradStepSize: The displacement step size used for construction of the numerical gradient (default 0.05)
   *        - numHessStepSize: The displacement step size used for construction of the numerical hessian (default 0.05)
   */
  HessianTaskSettings settings;

 private:
  const std::vector<std::shared_ptr<SystemController>> _activeSystems;
  const std::vector<std::shared_ptr<SystemController>> _passiveSystems;
};

} /* namespace Serenity */

#endif /* TASKS_HESSIANTASK_H_ */
