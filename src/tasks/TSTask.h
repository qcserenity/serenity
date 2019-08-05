/**
 * @file TSTask.h
 *
 * @date May 22, 2017
 * @author Jan Unsleber
 * @copyright \n
 * This file is part of the program Serenity.\n\n
 * Serenity is free software: you can redistribute it and/or modify
 * it under the terms of the LGNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.\n\n
 * Serenity is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.\n\n
 * You should have received a copy of the LGNU Lesser General
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
  TSTaskSettings():
    lst(false),
    lstqstonly(false),
    normalmode(1){
  };
  REFLECTABLE(
    (bool) lst,
    (bool) lstqstonly,
    (int) normalmode
  )
};
/**
 * @class  TSTask TSTask.h
 * @brief  A task for transition-state geometry-optimizations.
 */
class TSTask : public Task {
public:
  /**
   * @param systemController which provides you with all needed information about your system and
   *                         configuration
   */
  TSTask(std::shared_ptr<SystemController> ts,
            const std::vector<std::shared_ptr<SystemController> > env);
  virtual ~TSTask() = default;

  void run();

  ///@brief The settings for this task.
  TSTaskSettings settings;
private:
  std::shared_ptr<SystemController> _ts;
  std::vector<std::shared_ptr<SystemController> > _env;
};

} /* namespace QCpack */

#endif /* TSTASK_H_ */
