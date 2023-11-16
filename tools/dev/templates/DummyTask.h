/**
 * @file   DummyTask.h
 *
 * @date   Aug 12, 2014
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
#ifndef DUMMYTASK_H_
#define DUMMYTASK_H_
/* Includes */
#include "tasks/Task.h"
#include "settings/Reflection.h"

#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;
struct DummyTaskSettings {
  DummyTaskSettings(){

  };
  REFLECTABLE(
    (int) dummyoption
  )
};
/**
 * @class  DummyTask DummyTask.h
 * @brief  A placeholder for quickly making some functionality running temporarily. May be changed
 *         at any time for any purpose.
 */
class DummyTask : public Task {
public:
  /**
   * @param systemController A vector of all SystemControllers, which provide you with all needed
   *                         information about your systems and their settings.
   */
  DummyTask(const std::vector<std::shared_ptr<SystemController> > systemController);
  virtual ~DummyTask() = default;

  void run();
  /// @brief The vector of systemControllers to perform a dummy task on. Public-> anything may happen!
  std::vector<std::shared_ptr<SystemController> > _systemController;

  ///@brief The settings for this task.
  DummyTaskSettings settings;
};

} /* namespace QCpack */

#endif /* DUMMYTASK_H_ */
