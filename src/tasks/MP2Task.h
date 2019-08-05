/**
 * @file   MP2Task.h
 *
 * @date   Jul 14, 2014
 * @author Jan Unselber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef MP2TASK_H_
#define MP2TASK_H_
/* Include Serenity Internal Headers */
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
  MP2TaskSettings() :
   ri(true){
  };
  REFLECTABLE(
      (bool) ri
  )
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
   */
  MP2Task(std::shared_ptr<SystemController> systemController);

  /**
   * @brief Default destructor.
   */
  virtual ~MP2Task() = default;
  /**
   * @see Task
   */
  void run();

  /**
   * @brief The settings/keywords for the MP2Task: \n
   *        - ri: Use the RI approximation (default true).
   */
  MP2TaskSettings settings;
private:
  /// @brief The system controller.
  std::shared_ptr<SystemController> _systemController;

};

} /* namespace Serenity */

#endif /* MP2TASK_H_ */
