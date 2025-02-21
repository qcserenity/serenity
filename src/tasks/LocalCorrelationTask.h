/**
 * @file LocalCorrelationTask.h
 *
 * @date Apr 4, 2021
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

#ifndef SERENITY_LOCALCORRELATIONTASK_H
#define SERENITY_LOCALCORRELATIONTASK_H

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/CorrelatedMethodsOptions.h"
#include "settings/Reflection.h"
#include "tasks/LocalizationTask.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;
struct LocalCorrelationTaskSettings {
  LocalCorrelationTaskSettings() : maxCycles(100), normThreshold(1e-5), writePairEnergies(false) {
    lcSettings.enforceHFFockian = true;
    loc.splitValenceAndCore = true;
    lcSettings.linearScalingSigmaVector = true;
  }
  REFLECTABLE((unsigned int)maxCycles, (double)normThreshold, (bool)writePairEnergies)
 public:
  LocalCorrelationSettings lcSettings;
  LocalizationTaskSettings loc;
};
/**
 * @class
 * @brief A task that perform a local correlation calculation. It serves solely as
 *        an input simplification by calling the appropriate other tasks.
 */
class LocalCorrelationTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system              The system controller.
   * @param environmentSystems  The environment systems.
   *
   * Assigns private variables only.
   */
  LocalCorrelationTask(std::shared_ptr<SystemController> system,
                       std::vector<std::shared_ptr<SystemController>> environmentSystems = {});
  /**
   * @brief Default destructor
   */
  ~LocalCorrelationTask();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(LocalCorrelationTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.lcSettings.visitAsBlockSettings(v, blockname))
      return;
    if (c.loc.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in LocalCorrelationTaskSettings: " + blockname);
  }
  /**
   * @brief The task settings.
   * - maxCycles           Maximum number of cycles for the resiudal equations.
   * - normThreshold       Convergence threshold for the residual equations (energy and abs. max. residual)
   * - writePairEnergies   (Coupled cluster only) Write the pair energies to file.
   * - lcSettings          The local correlation settings.
   * - loc                 The localization task settings.
   */
  LocalCorrelationTaskSettings settings;
  /**
   * @see Task
   */
  void run();

  /**
   * @brief Getter for the correlation energy.
   * @param activeSystemController The system controller.
   * @param pnoMethod The local correlation method to get the energy for.
   * @return The energy.
   */
  static double getCorrelationEnergy(std::shared_ptr<SystemController> activeSystemController, Options::PNO_METHOD pnoMethod);

 private:
  // The system controller.
  std::shared_ptr<SystemController> _system;
  // Optional environment systems to be used in the Fock operator construction.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
};

} /* namespace Serenity */
#endif // SERENITY_LOCALCORRELATIONTASK_H
