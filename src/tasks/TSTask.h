/**
 * @file TSTask.h
 *
 * @date May 22, 2017
 * @author Jan Unsleber
 * @copyright \n
 * This file is part of the program Serenity.\n\n
 * Serenity is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.\n\n
 * Serenity is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.\n\n
 * You should have received a copy of the GNU Lesser General
 * Public License along with Serenity.
 * If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef TSTASK_H_
#define TSTASK_H_
/* Include Serenity Internal Headers */
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;
struct TSTaskSettings {
  TSTaskSettings() : lst(false), lstqstonly(false), normalmode(1){};
  REFLECTABLE((bool)lst, (bool)lstqstonly, (int)normalmode)
};
/**
 * @class  TSTask TSTask.h
 * @brief  A task for transition-state geometry-optimizations.
 */
class TSTask : public Task {
 public:
  /**
   * @param ts a shared pointer on the SystemController, which provides you with all needed information about the system
   * to be optimized
   * @param env a constant vector of shared pointers onto the SystemController, which provides you with all needed
   * information about the systems that represent both minima of the potential-energy surface
   */
  TSTask(std::shared_ptr<SystemController> ts, const std::vector<std::shared_ptr<SystemController>> env);
  virtual ~TSTask() = default;

  /**
   * @brief the function that actually runs the TSTask
   */
  void run();

  /**
   * @brief The settings/keywords for the TSTask:
   *        - DEFAULT (given one guess structure): does Bofill TS-search with given guess structure
   *        - DEFAULT (given guess TS and structures of minima): does QST and afterwards Bofill TS-search
   *
   *        - lst: specifies if LST is performed before QST
   *        - lstqstonly: specifies that only LST/QST should be done and NO Bofill TS-search afterwards
   *        - normalmode: specifies the reaction coordinate along which the TS is searched
   */

  TSTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _ts;
  std::vector<std::shared_ptr<SystemController>> _env;
};

} /* namespace Serenity */

#endif /* TSTASK_H_ */
