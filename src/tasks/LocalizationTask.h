/**
 * @file   LocalizationTask.h
 *
 * @date   Apr 22, 2014
 * @author Thomas Dresselhaus
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
#ifndef LOCALIZATIONTASK_H_
#define LOCALIZATIONTASK_H_
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
struct LocalizationTaskSettings {
  LocalizationTaskSettings():
    locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
    maxSweeps(1000){
  };
  REFLECTABLE(
    (Options::ORBITAL_LOCALIZATION_ALGORITHMS) locType,
    (unsigned int) maxSweeps
  )
};
/**
 * @class  LocalizationTask LocalizationTask.h
 * @brief  Localize orbitals.
 *
 * In general that means: perform unitary transformations among the occupied orbitals
 * to optimize some (localization) criterion.
 */
class LocalizationTask : public Task {
public:
  /**
   * @param system The system from which orbitals to localize are taken.
   */
  LocalizationTask(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Default destructor.
   */
  virtual ~LocalizationTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief The settings/keywords for LocalizationTask:
   *        - locType: The localization algorithm. The following algorihtm can be chosen:
   *          - IBO : Intrinsic Bond Orbitals (default)
   *          - PM : Pipek-Mezey
   *          - FB : Foster-Boys
   *          - ER : Edminston-Ruedenberg
   *        - maxSweeps : Maximum number of micro-iterations (default: 1000)
   */
  LocalizationTaskSettings settings;

private:
  template<Options::SCF_MODES T> void runByLastSCFMode();
  
  const std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* LOCALIZATIONTASK_H_ */
