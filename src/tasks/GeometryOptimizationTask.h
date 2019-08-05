/**
 * @file   GeometryOptimizationTask.h
 *
 * @date   May 21, 2015
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
#ifndef CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_
#define CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>


namespace Serenity {
/* Forward declarations */
class SystemController;
class OptimizerFactory;
class Optimizer;
using namespace Serenity::Reflection;
struct GeometryOptimizationTaskSettings {
  GeometryOptimizationTaskSettings():
    gradType(Options::GRADIENT_TYPES::ANALYTICAL),
    maxCycles(100),
    rmsgradThresh(1.0e-4),
    energyChangeThresh(5.0e-6),
    maxGradThresh(3.0e-4),
    stepThresh(2.0e-3),
    maxStepThresh(4.0e-3),
    numGradStepSize(0.001),
    printLevel(3),
    transInvar(false),
    naddKinFunc(Options::KINFUNCTIONALS::PW91K),
    naddXCFunc(Options::XCFUNCTIONALS::PW91),
    dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE),
    FaTmaxCycles(50),
    FaTConvThresh(1.0e-6),
    FaTgridCutOff(-1.0){}
  REFLECTABLE(
    (Options::GRADIENT_TYPES) gradType,
    (unsigned int) maxCycles,
    (double) rmsgradThresh,
    (double) energyChangeThresh,
    (double) maxGradThresh,
    (double) stepThresh,
    (double) maxStepThresh,
    (double) numGradStepSize,
    (unsigned int) printLevel,
    (bool) transInvar,
    (Options::KINFUNCTIONALS) naddKinFunc,
    (Options::XCFUNCTIONALS) naddXCFunc,
    (Options::DFT_DISPERSION_CORRECTIONS) dispersion,
    (unsigned int) FaTmaxCycles,
    (double) FaTConvThresh,
    (double) FaTgridCutOff
  )
};
/**
 * @class  GeometryOptimizationTask GeometryOptimizationTask.h
 * @brief  Task to optimize the geometry of the system
 */
template<Options::SCF_MODES SCFMode>
class GeometryOptimizationTask : public Task {
public:
  /**
   * @param systemController
   */

  GeometryOptimizationTask(
	      const std::vector<std::shared_ptr<SystemController> >& activeSystems,
	      const std::vector<std::shared_ptr<SystemController> >& passiveSystems
			= std::vector<std::shared_ptr<SystemController> >());
  /**
   * @brief Default destructor.
   */
  virtual ~GeometryOptimizationTask() = default;
  /**
   * @see Task
   */
  void run();

  /**
   * @brief The settings/keywords for the GeometryOptimizationTask:
   *        - gradType: The type of the gradient. The following options are available:
   *          - NUMERICAL
   *          - ANALYTICAL (default)
   *        - maxCycles: Maximum number of steps (default: 100)
   *        - gradNormThresh: The gradient convergence threshold (default: 1.0e-5)
   *        - numGradStepSize: The displacement step size used for construction of the numerical gradients (default 0.05)
   *        - printLevel:
   *          - 0 : Print nothing (default)
   *          - 1 : Print SCF results
   *          - 2 : Print SCF results and SCF cycle information
   *        - transInvar Make gradients translationally invariant.
   *        - naddKinFunc The non-additive kinetic energy functional.
   *        - naddXCFunc The non-addtive exchange--correlation functional.
   *        - dispersion The intersubsystem dispersion interaction.
   *        - FaTmaxCycles The maximum number of FaT iterations.
   *        - FaTConvThresh Convergence threshold for the FaT calculation.
   *        - FaTgridCutOff A distance cutoff for the integration grid used to calculate
   *                        the non-additive energy functional potentials.
   *
   */
  GeometryOptimizationTaskSettings settings;
private:
  std::vector<std::shared_ptr<SystemController> > _activeSystems;
  std::vector<std::shared_ptr<SystemController> > _passiveSystems;
};

} /* namespace Serenity */

#endif /* CONFIGURATION_TASKS_GEOMETRYOPTIMIZATIONTASK_H_ */
