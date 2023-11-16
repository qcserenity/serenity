/**
 * @file ExportGridTask.h
 *
 * @date: Mar 24, 2016
 * @author: Jan Unsleber
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

#ifndef TASKS_EXPORTGRIDTASK_H_
#define TASKS_EXPORTGRIDTASK_H_
/* Include Serenity Internal Headers */
#include "settings/Reflection.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace ::Serenity::Reflection;
struct ExportGridTaskSettings {
  ExportGridTaskSettings() : withAtomInfo(false){};
  REFLECTABLE((bool)withAtomInfo)
};
/**
 * @class ExportGridTask ExportGridTask.h
 * @brief A Task to export grids in ASCII format.
 */
class ExportGridTask : public Task {
 public:
  /**
   * @brief Constructor.
   *
   * @param systemController The system.
   */
  ExportGridTask(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Default destructor.
   */
  virtual ~ExportGridTask() = default;

  /**
   * @see Task
   */
  void run();

  /**
   * @brief The settings/keywords for ExportGridTask:
   *        -withAtomInfo: If true, additionally print information about the atoms (default: false),
   *            please first check if gridPointSorting is set to false in the systems' grid block
   */
  ExportGridTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* TASKS_EXPORTGRIDTASK_H_ */
