/**
 * @file SystemSplittingTask.h
 *
 * @author Moritz Bensberg
 * @date Jan 8, 2020
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

#ifndef TASKS_SYSTEMSPLITTINGTASK_H_
#define TASKS_SYSTEMSPLITTINGTASK_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h" //Output praefixes.
#include "settings/MiscOptions.h"   //SYSTEM_SPLITTING_ALGORITHM.
#include "tasks/Task.h"             //inherits from.
/* Include Std and External Headers */
#include <Eigen/Dense> //Eigen::VectorXi
#include <memory>      //smrt_ptr.
#include <string>      //Output praefixes.
#include <vector>      //std::vector.

namespace Serenity {

/* Forward Declarations */
class SystemController;

struct SystemSplittingTaskSettings {
  SystemSplittingTaskSettings()
    : orbitalThreshold(0.4), systemPartitioning(Options::SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH) {
  }
  REFLECTABLE((double)orbitalThreshold, (Options::SYSTEM_SPLITTING_ALGORITHM)systemPartitioning)
 public:
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("SPLIT")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};

/**
 * @class SystemSplittingTask SystemSplittingTask.h
 * @brief A task that partitions a supersystem into subsystems.
 */
template<Options::SCF_MODES SCFMode>
class SystemSplittingTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param supersystem The supersystem to be partitioned.
   * @param subsystems The subsystems to be filled.
   */
  SystemSplittingTask(std::shared_ptr<SystemController> supersystem, std::vector<std::shared_ptr<SystemController>> subsystems);
  virtual ~SystemSplittingTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief The task settings.
   *   systemPartitioning: The partitioning algorithm.
   *   orbitalThreshold:   Threshold for population-threshold orbital selection for
   *                       subsystem 0.
   */
  SystemSplittingTaskSettings settings;
  /**
   * @brief Getter for the final orbital partitioning.
   *   The assignment vector contains the subsystem index for each occupied orbital.
   * @return The assignment vector.
   */
  const SpinPolarizedData<SCFMode, Eigen::VectorXi>& getFinalAssignment();

 private:
  // The supersystem.
  std::shared_ptr<SystemController> _supersystem;
  // The subsystems.
  std::vector<std::shared_ptr<SystemController>> _subsystems;
  // Helper function to search for the subsystem atoms in the supersystem geometry.
  std::vector<unsigned int> findAtoms(std::shared_ptr<SystemController> subsystem);
  // The final assignments. Only available after executing run().
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXi>> _assignment = nullptr;
  // Helper function to make sure that the supersystem can be partitioned into the given fragments.
  void checkInput();
};

} /* namespace Serenity */

#endif /* TASKS_SYSTEMSPLITTINGTASK_H_ */
