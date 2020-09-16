/**
 * @file   MP2Task.h
 *
 * @date   Jul 14, 2014
 * @author Jan Unselber
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
#ifndef MP2TASK_H_
#define MP2TASK_H_
/* Include Serenity Internal Headers */
#include "memory/MemoryManager.h"                               //Memory handling
#include "misc/SerenityError.h"                                 //Error messages.
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //LocalCorrelationSettings.
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;

struct MP2TaskSettings {
  MP2TaskSettings() : mp2Type(Options::MP2_TYPES::RI), maxResidual(1e-5), maxCycles(100) {
    lcSettings.pnoSettings = Options::PNO_SETTINGS::NORMAL;
  };
  REFLECTABLE((Options::MP2_TYPES)mp2Type, (double)maxResidual, (unsigned int)maxCycles)
 public:
  LocalCorrelationSettings lcSettings;
};
/**
 * @class MP2Task MP2Task.h
 * @brief Performs an MP2 energy correction with the given orbitals.
 */
template<Options::SCF_MODES SCFMode>
class MP2Task : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController
   * @param environmentSystems For local MP2 the occupied env. orbitals are projected out of
   *                           the space of the virtual orbitals.
   */
  MP2Task(std::shared_ptr<SystemController> systemController, std::vector<std::shared_ptr<SystemController>> environmentSystems);

  /**
   * @brief Default destructor.
   */
  virtual ~MP2Task() = default;
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
  void visit(MP2TaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
    }
    else if (!c.lcSettings.visitSettings(v, blockname)) {
      throw SerenityError((string) "Unknown settings block in CoupledClusterTaskSettings: " + blockname);
    }
  }

  /**
   * @brief The settings/keywords for the MP2Task: \n
   *        - mp2Type [RI]       : MP2-Type.
   *        - maxResidual [1e-5] : (Local only) Maximum tolerated residual.
   *        - maxCycles [100]    : (Local only) Maximum number of amplitude optimization cycles.
   *        -lcSettings Local correlation settings. @see LocalCorrelationController.h
   *
   */
  MP2TaskSettings settings;

 private:
  /// @brief The system controller.
  std::shared_ptr<SystemController> _systemController;
  /// @brief The environment systems.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
};

} /* namespace Serenity */

#endif /* MP2TASK_H_ */
