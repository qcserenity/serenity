/**
 *
 * @file   DispersionCorrectionTask.h
 *
 * @date   Dec 02, 2015
 * @author Jan Unsleber
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

#ifndef DISPERSIONCORRECTIONTASK_H_
#define DISPERSIONCORRECTIONTASK_H_
/* Include Serenity Internal Headers */
#include "settings/DFTOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class SystemController;

using namespace Serenity::Reflection;
struct DispersionCorrectionTaskSettings {
  DispersionCorrectionTaskSettings()
    : dispType(Options::DFT_DISPERSION_CORRECTIONS::D3BJ),
      functional(CompositeFunctionals::XCFUNCTIONALS::NONE),
      gradient(false),
      hessian(false){};
  REFLECTABLE((Options::DFT_DISPERSION_CORRECTIONS)dispType, (CompositeFunctionals::XCFUNCTIONALS)functional,
              (bool)gradient, (bool)hessian)
};
/**
 * @class DispersionCorrectionTask DispersionCorrectionTask.h
 * @brief Calculates and displays the DFT dispersion correction for a given system.
 */
class DispersionCorrectionTask : public Task {
 public:
  /**
   * @brief Constructor
   * @param systemConstroller The system for which the correction shall be calculated.
   */
  DispersionCorrectionTask(std::shared_ptr<SystemController> systemController);

  /**
   * @brief Default destructor
   */
  virtual ~DispersionCorrectionTask() = default;

  void run() override final;

  /**
   * @brief The settings/keywords for DispersionCorrectionTask: \n
   *        - dispType : The type of dispersion correction to be used. The following corrections are available:
   *           - none     : No correction
   *           - D3       : The third set of parameters (with zero damping, also called D3(0)).
   *           - D3ABC    : The third set of parameters (with zero damping, also called D3(0)) and 3 center correction
   * term.
   *           - D3BJ     : The third set of parameters with Becke-Johnson damping (default).
   *           - D3BJiABC : The third set of parameters with Becke-Johnson damping and 3 center correction term.
   *        - functional : The used XC-functional (Default: The functional which is defined in the system).
   *        - gradient: If true, calculate dispersion gradient contribution (default: false).
   *        - hessian: If true, calculate dispersion hessian contribution (default: false).
   */
  DispersionCorrectionTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* DISPERSIONCORRECTIONTASK_H_ */
