/**
 * @file CoupledClusterTask.h
 *
 * @date Apr 3, 2016
 * @author Jan Unsleber
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

#ifndef TASKS_COUPLEDCLUSTERTASK_H_
#define TASKS_COUPLEDCLUSTERTASK_H_

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

struct CoupledClusterTaskSettings {
  CoupledClusterTaskSettings():
    level(Options::CC_LEVEL::CCSD),
	maxCycles(100),
	normThreshold(1e-7){}
  REFLECTABLE(
    (Options::CC_LEVEL) level,
	(unsigned int) maxCycles,
	(double) normThreshold
  )
};



/**
 * @class CoupledClusterTask CoupledClusterTask.h
 * @brief A class to run coupled cluster calculations.
 */
class CoupledClusterTask : public Task {
public:
  /**
   * @brief Constructor
   * @param systemController The system of interest.
   */
  CoupledClusterTask(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Destructor.
   */
  virtual ~CoupledClusterTask() = default;

  /**
   * @brief @see Task
   */
  void run() override final;

  /**
  * @brief Settings/Keywords for CoupledClusterTask: \n
  *        -CC_LEVEL: The excitation level used in the CC calculation. The following methods are implemented:
  *          - CCSD (default)
  *          - CCSD_T
  */
  CoupledClusterTaskSettings settings;

private:
  std::shared_ptr<SystemController> _systemController;
};



} /* namespace Serenity */

#endif /* TASKS_COUPLEDCLUSTERTASK_H_ */
