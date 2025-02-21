/**
 * @file   PopAnalysisTask.h
 *
 * @date   last rework Jun 30, 2017
 * @author Thomas Dresselhaus, last rework Jan Unsleber
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
#ifndef POPANALYSISTASK_H_
#define POPANALYSISTASK_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/LocalizationOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;

struct PopAnalysisTaskSettings {
  PopAnalysisTaskSettings() : algorithm(Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN) {
  }
  REFLECTABLE((Options::POPULATION_ANALYSIS_ALGORITHMS)algorithm)
};

/**
 * @class  PopAnalysisTask PopAnalysisTask.h
 * @brief  Task to perform a population analysis. Currently fixed to Mulliken.
 */
template<Options::SCF_MODES SCFMode>
class PopAnalysisTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param systemController from which e.g. the active density is taken to calculate the populations/charges.
   */
  PopAnalysisTask(std::shared_ptr<SystemController> systemController);
  /**
   * @brief Default destructor.
   */
  virtual ~PopAnalysisTask() = default;

  void run();

  /**
   * @brief The settings/keywords for PopAnalysisTask:
   *          -Mulliken : Perform Mulliken population analysis (default)
   *          -Hirschfeld : Perform Hirschfeld population analysis
   *          -Becke : Perform Becke population analysis
   *          -IAO : Perform IAO population analysis.
   *          -IAOShell : Perform Shell-wise IAO population analysis.
   *          -CM5 : Perform CM5 population analysis.
   *          -CHELPG: Perform a CHELPG population analysis.
   */
  PopAnalysisTaskSettings settings;

 private:
  void print(std::string type, const SpinPolarizedData<SCFMode, Eigen::VectorXd>& populations);

  const std::shared_ptr<SystemController> _systemController;
  std::string _type;
  SpinPolarizedData<SCFMode, std::string> _modestring;
};

} /* namespace Serenity */

#endif /* POPANALYSISTASK_H_ */
