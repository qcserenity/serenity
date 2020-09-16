/**
 * @file EvaluateEnergyTask.h
 *
 * @author Moritz Bensberg
 * @date Mar 10, 2020
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

#ifndef TASKS_EVALUATEENERGYTASK_H_
#define TASKS_EVALUATEENERGYTASK_H_

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Local correlation settings for local MP2.
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward Declarations */
class SystemController;

using namespace Serenity::Reflection;
struct EvaluateEnergyTaskSettings {
  EvaluateEnergyTaskSettings() : mp2Type(Options::MP2_TYPES::RI), maxResidual(1e-5), maxCycles(100) {
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = PNO_METHOD::DLPNO_MP2;
  }
  REFLECTABLE((Options::MP2_TYPES)mp2Type, (double)maxResidual, (int)maxCycles)
 public:
  LocalCorrelationSettings lcSettings;
};
/**
 * @class EvaluateEnergyTask EvaluateEnergyTask.h
 * @brief Evaluates the energy for the given system, density and settings.
 */
template<Options::SCF_MODES SCFMode>
class EvaluateEnergyTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController The system for which the energy should be evaluated.
   */
  EvaluateEnergyTask(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Default destructor.
   */
  virtual ~EvaluateEnergyTask() = default;

  /**
   * @see Task
   */
  void run();
  /**
   * brief Settings.
   *        -mp2Type:               Type of MP2 used for double-hybrid correlation part.
   *        -maxResidual:           Maximum residual threshold for local MP2.
   *        -maxCycles:             Maximum number of iterations before cancelling the amplitude optimization
   *                                in local MP2.
   *        -lcSettings:            Local correlation settings for local MP2.
   */
  EvaluateEnergyTaskSettings settings;
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(EvaluateEnergyTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
    }
    else if (!c.lcSettings.visitSettings(v, blockname)) {
      throw SerenityError((string) "Unknown settings block in EvaluateEnergyTaskSettings: " + blockname);
    }
  }

 private:
  // The system controller.
  std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* TASKS_EVALUATEENERGYTASK_H_ */
