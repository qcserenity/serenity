/**
 * @file   Task.h
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef TASK_H_
#define TASK_H_
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {

extern Options::GLOBAL_PRINT_LEVELS GLOBAL_PRINT_LEVEL;
extern IOOptions iOOptions;

class SystemController;

using namespace Serenity::Reflection;
struct GeneralTaskSettings {
  GeneralTaskSettings() : printLevel(Options::GLOBAL_PRINT_LEVELS::NORMAL), timingsPrintLevel(1){};
  REFLECTABLE((Options::GLOBAL_PRINT_LEVELS)printLevel, (unsigned int)timingsPrintLevel)
};

/**
 * @class Task Task.h
 * @brief Interface for jobs that can actually be performed.
 *
 * A Task can be set up in any possible way and is, if part of the task list, triggered during the
 * program run.
 */
class Task {
 public:
  Task() = default;
  virtual ~Task() = default;
  /**
   * This method invokes the actual task. Everything may happen here.
   */
  virtual void run() = 0;

  /**
   * @brief Parse general settings to global settings for the output.
   *
   *        This function is called in the main (serenity.cpp) directly
   *        before executing the run() function for the task.
   */
  void parseGeneralSettings() {
    GLOBAL_PRINT_LEVEL = generalSettings.printLevel;
    iOOptions.timingsPrintLevel = generalSettings.timingsPrintLevel;
  }

  /**
   * @brief Settings in common by every task.
   *
   *   printLevel:         Print level used during task execution.
   *   timingsPrintLevel:  Timings print level used during task execution.
   *
   */

  GeneralTaskSettings generalSettings;

  /**
   * @brief Parses the settings to the task. Can be over defined by a new implementation
   *        in the respective task.
   * @param c The settings of any kind.
   * @param v The visitor.
   * @param blockname The name of the input block.
   */
  template<class C, class Visitor>
  void visit(C& c, Visitor v, std::string blockname) {
    (void)blockname;
    visit_each(c, v);
  }

  /**
   * @brief Determines which SCFMode a task should be run in and changes the SCFMode of all active systems and
   * supersystems to this new SCFMode.
   * @return The new SCFMode.
   */
  void avoidMixedSCFModes(Options::SCF_MODES runMode, std::vector<std::shared_ptr<SystemController>> activeSystems,
                          std::vector<std::shared_ptr<SystemController>> superSystems);

  /**
   * @brief Determines which SCFMode a task should be run in and changes the SCFMode of all active systems and
   * supersystems to this new SCFMode.
   * @return The new SCFMode.
   */
  void avoidMixedSCFModes(Options::SCF_MODES runMode, std::vector<std::shared_ptr<SystemController>> activeSystems);

  /**
   * @brief Determines which SCFMode a task should be run in and changes the SCFMode of all active systems and
   * supersystems to this new SCFMode.
   * @return The new SCFMode.
   */
  void avoidMixedSCFModes(Options::SCF_MODES runMode, std::shared_ptr<SystemController> activeSystem);
};

} /* namespace Serenity */

#endif /* TASK_H_ */
