/**
 * @file   PCMInteractionEnergyTask.h
 *
 * @date   Nov 12, 2020
 * @author Moritz Bensberg
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
#ifndef PCMINTERACTIONENERGYTASK_H_
#define PCMINTERACTIONENERGYTASK_H_
/* Include Serenity Internal Headers */
#include "settings/PCMSettings.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;
struct PCMInteractionEnergyTaskSettings {
  PCMInteractionEnergyTaskSettings() : dummy(0) {
    pcm.use = true;
  }
  REFLECTABLE((int)dummy // We need at least one setting.
  )
 public:
  PCMSettings pcm;
};

/**
 * @class  PCMInteractionEnergyTask PCMInteractionEnergyTask.h
 * @brief A class that allows the calculation of the dielectric energy correction from
 *        PCM for a given system (with associated density). Note that the apparent surface
 *        charges are not updated self-consistently!
 *
 *        The PCMSettings used may be specified via the task input.
 *
 *        TODO: Allow the user to specify environment systems in order to modify the cavity
 *              shape.
 *        TODO: Manual.
 */
class PCMInteractionEnergyTask : public Task {
 public:
  /**
   * @param systemController The system.
   */
  PCMInteractionEnergyTask(std::shared_ptr<SystemController> systemController);
  /**
   * |brief Default destructor
   */
  virtual ~PCMInteractionEnergyTask() = default;
  /*
   * @brief Execute the task.
   */
  void run();

  ///@brief The settings for this task.
  PCMInteractionEnergyTaskSettings settings;
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(PCMInteractionEnergyTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.pcm.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in PCMInteractionEnergyTaskSettings: " + blockname);
  }

 private:
  /// @brief The systemController
  std::shared_ptr<SystemController> _systemController;
  template<Options::SCF_MODES SCFMode>
  void runSpinPolarized();
};

} // namespace Serenity

#endif /* PCMINTERACTIONENERGYTASK_H_ */
