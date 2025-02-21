/**
 * @file   GradientTask.h
 *
 * @date   Mar 23, 2015
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
#ifndef CONFIGURATION_TASKS_GRADIENTTASK_H_
#define CONFIGURATION_TASKS_GRADIENTTASK_H_
/* Include Serenity Internal Headers */
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
#include "potentials/bundles/PotentialBundle.h"
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
struct GradientTaskSettings {
  GradientTaskSettings()
    : gradType(Options::GRADIENT_TYPES::ANALYTICAL),
      numGradStepSize(0.001),
      transInvar(false),
      FaTmaxCycles(50),
      FaTenergyConvThresh(1.0e-6),
      FDEgridCutOff(-1.0),
      print(true),
      printTotal(false) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::PW91K;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
  }
  REFLECTABLE((Options::GRADIENT_TYPES)gradType, (double)numGradStepSize, (bool)transInvar, (int)FaTmaxCycles,
              (double)FaTenergyConvThresh, (double)FDEgridCutOff, (bool)print, (bool)printTotal)
 public:
  EmbeddingSettings embedding;
};
/**
 * @class  GradientTask GradientTask.h
 * @brief  Task to calculate cartesian gradients of the system
 */
template<Options::SCF_MODES SCFMode>
class GradientTask : public Task {
 public:
  /**
   * @param systemController         The system for which gradients shall be calculated
   */
  GradientTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
               const std::vector<std::shared_ptr<SystemController>>& passiveSystems =
                   std::vector<std::shared_ptr<SystemController>>());
  /**
   * @brief Default destructor.
   */
  virtual ~GradientTask();
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
  void visit(GradientTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    // If reached, the blockname is unknown.
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    throw SerenityError((std::string) "Unknown block in FreezeAndThawTaskSettings: " + blockname);
  }
  /**
   * @brief The settings/keywords for GradientTask:
   *        - gradType: The type of the gradient. The following options are available:
   *          - NUMERICAL
   *          - ANALYTICAL (default)
   *        - numGradStepSize: The displacement step size used for construction of the numerical gradients (default 0.05)
   */
  GradientTaskSettings settings;

 private:
  const std::vector<std::shared_ptr<SystemController>> _activeSystems;
  const std::vector<std::shared_ptr<SystemController>> _passiveSystems;
  void printTotalGradient();
  void printPointChargeGradients();
};

} /* namespace Serenity */

#endif /* CONFIGURATION_TASKS_GRADIENTTASK_H_ */
