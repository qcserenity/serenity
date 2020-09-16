/**
 *  @file   EDATask.h
 *
 *  @date   Mar 29, 2016
 *  @author Moritz Bensberg
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

#ifndef TASKS_EDATASK_H_
#define TASKS_EDATASK_H_

/* Include Serenity Internal Headers */
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declaration */
class SystemController;
using namespace Serenity::Reflection;
namespace Options {
enum class SCF_MODES;
}

struct EDATaskSettings {
  EDATaskSettings() : dummy(false){};
  REFLECTABLE((bool)dummy)
};

/**
 * @class EDATask EDATask.h
 * @brief The task for the calculation and printing of a energy decomposition scheme for molecular
 *        interaction within the Hartree-Fock Approximation, according as
 *   Ref: Kitaura-Morokuma : K. Kitaura and K. Morokuma, Int. J. Quantum Chem. 10, 325 (1976)
 */
template<Options::SCF_MODES SCFMode>
class EDATask : public Task {
 public:
  /**
   * @brief The constructor.
   * @param systems A list of the systems to be used. (Only the first tow are used.)
   */
  EDATask(std::vector<std::shared_ptr<SystemController>> systems);

  /**
   * @brief the "run" routine which calls the calculator.
   */
  void run();

  /**
   * @brief The settings/keywords for the EDATask.\n
   */
  EDATaskSettings settings;

  /**
   * @brief default destructor
   */
  virtual ~EDATask() = default;

 private:
  /* the systems */
  /// @brief First system.
  std::shared_ptr<SystemController> _systemA;
  /// @brief Second system.
  std::shared_ptr<SystemController> _systemB;
};
} /* namespace Serenity */

#endif /* TASKS_EDATASK_H_ */
