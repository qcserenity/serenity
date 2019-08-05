/**
 * @file   ScfTask.h
 *
 * @date   Mar 7, 2014
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
#ifndef SCFTASK_H_
#define SCFTASK_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>


namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>class ElectronicStructureCalculator;
class SystemController;

using namespace Serenity::Reflection;

struct ScfTaskSettings {
  ScfTaskSettings():
    restart(false),
    fractionalDegeneracy(false){}
  REFLECTABLE(
      (bool) restart,
      (bool) fractionalDegeneracy
  )
};


/**
 * @class ScfTask ScfTask.h
 * @brief Performs an SCF calculation, i.e. an electronic structure calculation.
 */
template <Options::SCF_MODES SCFMode>class ScfTask : public Task {
public:
  /**
   * @param systemController
   */
  ScfTask(std::shared_ptr<SystemController> systemController);

  virtual ~ScfTask() = default;

  void run();

  /**
   * @brief The settings/keywords for SCFTask:
   *        -restart : Uses old orbitals to restart SCF (tries to restart from h5 file if available)
   */
  ScfTaskSettings settings;

  /**
   * @brief Getter for the electronic structure.
   * @return Returns the electronic structure.
   */
  std::shared_ptr<ElectronicStructureCalculator<SCFMode> > getElectronicStructureCalculator();

private:
  const std::shared_ptr<SystemController> _systemController;
};

} /* namespace Serenity */

#endif /* SCFTASK_H_ */
