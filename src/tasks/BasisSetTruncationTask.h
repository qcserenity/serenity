/**
 * @file BasisSetTruncationTask.h
 *
 * @date Sep 24, 2018
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

#ifndef TASKS_BASISSETTRUNCATIONTASK_H_
#define TASKS_BASISSETTRUNCATIONTASK_H_

/* Include Serenity Internal Headers */
#include "settings/MiscOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
using namespace Serenity::Reflection;
/* Forward declarations */
class SystemController;
class Atom;
class Geometry;
template<Options::SCF_MODES SCFMode>
class ElectronicStructure;

struct BasisSetTruncationTaskSettings {
  BasisSetTruncationTaskSettings()
    : truncAlgorithm(Options::BASIS_SET_TRUNCATION_ALGORITHMS::NET_POPULATION), netThreshold(1.0e-4), truncationFactor(0.0) {
  }
  REFLECTABLE((Options::BASIS_SET_TRUNCATION_ALGORITHMS)truncAlgorithm, (double)netThreshold, (double)truncationFactor)
 public:
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("TRUNC")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};

/**
 * @class BasisSetTruncationTask BasisSetTruncationTask.h
 * @brief Performs a truncation of the basis set of the given system. Only basis functions centered on dummy atoms
 *        are truncated. All others are considered to be the core basis of the system. Note that this task can
 *        manipulate the geometry of the given system such that no dummy atoms without any basis functions survive.
 *        This effectively truncates any associated fitting basis.
 */
template<Options::SCF_MODES SCFMode>
class BasisSetTruncationTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller.
   */
  BasisSetTruncationTask(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor.
   */
  virtual ~BasisSetTruncationTask() = default;
  /**
   * @brief Default run method.
   */
  void run();

  /**
   * -- truncAlgorithm            The truncation algorithm.
   * -- netThreshold              The threshold for the net-population truncation
   * -- truncationFactor          The factor for the primitive net-population truncation.
   */
  BasisSetTruncationTaskSettings settings;

 private:
  /// @brief The system controller.
  std::shared_ptr<SystemController> _system;
};

} /* namespace Serenity */

#endif /* TASKS_BASISSETTRUNCATIONTASK_H_ */
