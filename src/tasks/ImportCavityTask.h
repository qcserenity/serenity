/**
 * @file   ImportCavityTask.h
 *
 * @date   Nov 08, 2024
 * @author Lukas Paetow
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
#ifndef TASKS_IMPORTCAVITYTASK_H_
#define TASKS_IMPORTCAVITYTASK_H_
/* Include Serenity Internal Headers */
#include "geometry/Sphere.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
enum class MOLECULAR_SURFACE_TYPES;
using namespace Serenity::Reflection;
struct ImportCavityTaskSettings {
  ImportCavityTaskSettings() : cavityPath(""), vdwcavityPath(""), fdecavity(false){};
  REFLECTABLE((std::string)cavityPath, (std::string)vdwcavityPath, (bool)fdecavity)
};
/**
 * @class  ImportCavityTask ImportCavityTask.h
 * @brief  A task to load the data of a solvation model including the cavity and the optimized charges from a few HDF5
 * files.
 */
class ImportCavityTask : public Task {
 public:
  /**
   * @param systemController A SystemController, which provides you with all needed
   *                         information about your system and its settings.
   */
  ImportCavityTask(std::shared_ptr<SystemController> systemController);
  virtual ~ImportCavityTask() = default;

  SphereType stringToSphereType(std::string typestring);

  void run();

  // Handles the main work of this task.
  void importCavityFunction(std::string filepath, MOLECULAR_SURFACE_TYPES surfaceType);

  /// @brief The active system named in the input file for this task
  std::shared_ptr<SystemController> _systemController;

  /**
   * @brief The settings/keywords for ImportCavityTask:
   *        -cavityPath: The path of the molecular solvent cavity to be imported.
   *        -vdwcavityPath: The path of the molecular solvent Van der Waals cavity to be imported.
   *        -fdecavity: Set this to true if you are dealing with a cavity that surrounds several subsystems.
   */
  ImportCavityTaskSettings settings;
};

} // namespace Serenity

#endif /* TASKS_IMPORTCAVITYTASK_H_ */
