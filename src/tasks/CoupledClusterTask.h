/**
 * @file CoupledClusterTask.h
 *
 * @date Apr 3, 2016
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

#ifndef TASKS_COUPLEDCLUSTERTASK_H_
#define TASKS_COUPLEDCLUSTERTASK_H_

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/CorrelatedMethodsOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;

struct CoupledClusterTaskSettings {
  CoupledClusterTaskSettings()
    : level(Options::CC_LEVEL::CCSD), maxCycles(100), normThreshold(1e-5), writePairEnergies(false) {
    lcSettings.enforceHFFockian = true;
  }
  REFLECTABLE((Options::CC_LEVEL)level, (unsigned int)maxCycles, (double)normThreshold, (bool)writePairEnergies)
 public:
  LocalCorrelationSettings lcSettings;
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
  CoupledClusterTask(std::shared_ptr<SystemController> systemController,
                     std::vector<std::shared_ptr<SystemController>> environmentSystems = {});
  /**
   * @brief Destructor.
   */
  virtual ~CoupledClusterTask() = default;
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(CoupledClusterTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    // If reached, the blockname is unknown.
    if (c.lcSettings.visitAsBlockSettings(v, blockname))
      return;
    throw SerenityError((std::string) "Unknown settings block in CoupledClusterTaskSettings: " + blockname);
  }

  /**
   * @brief @see Task
   */
  void run() override final;

  /**
   * @brief Settings/Keywords for CoupledClusterTask: \n
   *        -CC_LEVEL: The excitation level used in the CC calculation. The following methods are implemented:
   *          - CCSD (default)
   *          - CCSD_T
   *          - DLPNO_CCSD
   *          - DLPNO_CCSD_T0
   *        -maxCycles
   *        -normThreshold
   *        -lcSettings Local correlation settings. @see LocalCorrelationController.h
   */
  CoupledClusterTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _systemController;

  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  // Performs canonical calculation.
  Eigen::VectorXd canonicalCalculation();
  // Performs local calculation.
  Eigen::VectorXd localCalculation();
};

} /* namespace Serenity */

#endif /* TASKS_COUPLEDCLUSTERTASK_H_ */
