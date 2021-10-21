/**
 * @file SystemAdditionTask.h
 *
 * @author Moritz Bensberg
 * @date Jan 9, 2020
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

#ifndef TASKS_SYSTEMADDITIONTASK_H_
#define TASKS_SYSTEMADDITIONTASK_H_

/* Include Serenity Internal Headers */
#include "tasks/Task.h" //inherits from.
/* Include Std and External Headers */
#include <memory> //smrt_ptr.
#include <vector> //std::vector.

namespace Serenity {

/* Forward Declarations */
class SystemController;
class Geometry;
struct SystemAdditionTaskSettings {
  SystemAdditionTaskSettings() : checkSuperGeom(false), checkSuperCharge(false), addOccupiedOrbitals(true) {
  }
  REFLECTABLE((bool)checkSuperGeom, (bool)checkSuperCharge, (bool)addOccupiedOrbitals)
 public:
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("ADD")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};
/**
 * @class SystemAdditionTask SystemAdditionTask.h
 * @brief A task that combines an arbitrary number of subsystems to a supersystem.
 */
template<Options::SCF_MODES SCFMode>
class SystemAdditionTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param supersystem The supersystem to be constructed.
   * @param subsystems The subsystems to be combined.
   */
  SystemAdditionTask(std::shared_ptr<SystemController> supersystem, std::vector<std::shared_ptr<SystemController>> subsystems);
  virtual ~SystemAdditionTask() = default;

  /**
   * @see Task
   */
  void run();
  /**
   * @brief The task settings.
   *   checkSuperGeom      Assert that the supersystem geometry is not changed.
   *   checkSuperCharge    Assert that spin and charge of the supersystem is not
   *                       changed.
   *   addOccupiedOrbitals Create an initial guess for the supersystem-electronic
   *                       structure from the occupied orbitals of the subsystems.
   */
  SystemAdditionTaskSettings settings;

 private:
  // The supersystem controller.
  std::shared_ptr<SystemController> _supersystem;
  // The subsystem controller.
  std::vector<std::shared_ptr<SystemController>> _subsystems;
  // Check the supersystem geometry for changes.
  void checkGeometry(std::shared_ptr<Geometry> supersystemGeometry);
};

} /* namespace Serenity */

#endif /* TASKS_SYSTEMADDITIONTASK_H_ */
