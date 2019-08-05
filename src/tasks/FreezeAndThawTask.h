/**
 * @file FreezeAndThawTask.h
 *
 * @date Oct 15, 2015
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

#ifndef CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_
#define CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "dft/Functional.h"
#include "data/OrbitalController.h"
#include "settings/Reflection.h"
#include "system/SystemController.h"
#include "settings/Options.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>


namespace Serenity {

using namespace Serenity::Reflection;

struct FreezeAndThawTaskSettings {
  FreezeAndThawTaskSettings():
    naddKinFunc(Options::KINFUNCTIONALS::PW91K),
    naddXCFunc(Options::XCFUNCTIONALS::PW91),
    embeddingMode(Options::KIN_EMBEDDING_MODES::NADD_FUNC),
    dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE),
    maxCycles(50),
    convThresh(1.0e-6),
    smoothFactor(0.0),
    potentialBasis(""),
    singValThreshold(0.0),
    lbDamping(0.995),
    lbCycles(0),
    carterCycles(0),
    gridCutOff(-1.0),
    printLevel(2),
    makeSuperSystemBasis(false),
    smallSupersystemGrid(false),
    basisExtThresh(5.0e-2),
    extendBasis(false),
    truncateProjector(false),
    projecTruncThresh(1.0e+1),
    distantKinFunc(false),
    useConvAcceleration(false),
    diisStart(5.0e-5),
    diisEnd(1e-4),
    levelShiftParameter(1e+6),
    basisFunctionRatio(0.0),
    borderAtomThreshold(0.02)
  {}
  REFLECTABLE(
      (Options::KINFUNCTIONALS) naddKinFunc,
      (Options::XCFUNCTIONALS) naddXCFunc,
      (Options::KIN_EMBEDDING_MODES) embeddingMode,
      (Options::DFT_DISPERSION_CORRECTIONS) dispersion,
      (unsigned int) maxCycles,
      (double) convThresh,
      (double) smoothFactor,
      (std::string) potentialBasis,
      (double) singValThreshold,
      (double) lbDamping,
      (unsigned int) lbCycles,
      (unsigned int) carterCycles,
      (double) gridCutOff,
      (unsigned int) printLevel,
      (bool) makeSuperSystemBasis,
      (bool) smallSupersystemGrid,
      (double) basisExtThresh,
      (bool) extendBasis,
      (bool) truncateProjector,
      (double) projecTruncThresh,
      (bool) distantKinFunc,
      (bool) useConvAcceleration,
      (double) diisStart,
      (double) diisEnd,
      (double) levelShiftParameter,
      (double) basisFunctionRatio,
      (double) borderAtomThreshold
  )
};

/**
 * @class FreezeAndThawTask FreezeAndThawTask
 * @brief A task for freeze and thaw cycles using only FDE (for now).
 */
template<Options::SCF_MODES SCFMode>
class FreezeAndThawTask : public Task {
public:
  /**
   * @brief Constructor.
   *
   * @param activeSystems          A list of all the active systems.
   * @param passiveSystems         A list of all the passive systems (never turning active).
   */
  FreezeAndThawTask(
      const std::vector<std::shared_ptr<SystemController> >& activeSystems,
      const std::vector<std::shared_ptr<SystemController> >& passiveSystems = {});
  /**
   * @brief Default destructor.
   */
  virtual ~FreezeAndThawTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief The settings/keywords for the FreezeAndThawTask:
   *        - naddKinFunc : The non-additive kinetic energy functional (default: TF)
   *        - naddXCFunc: The non-additive exchange-correlation functional (default: BP86)
   *        - maxCycles: The maximum number of FAT iterations (default: 5)
   *        - energyConvThresh : The convergence threshold for the energy (default: 1,0e-6)
   *        - gradConvThresh : The convergence threshold for the gradient (used for FAT geometry optimizations, default: 1.0e-5)
   *        - opt: If true, perform FAT geometry optimization (default: false)
   *        - exactNaddKin: If true, use reconstructed non-additive kinetic potential (ignores the naddKinFunc keyword, Default: false)
   *        - projectionOperator: Specifies which projection operator is used. (default: NONE)
   *        - gridCutOff :  Modifies the super-system grid used to evaluate the FDE potential. See FDETask.
   *        - basisExtensionThreshold: The threshold for the basis set extension (default: 5.0 e-2). See BasisExtension.h
   *        - extendBasis: A flag whether the basis should be extended (default: false). See BasisExtension.h
   *        - truncateProjector: A flag whether the projector should be truncated (default: false). See HuzinagaFDEProjectionPotential.h.
   *        - projectionTruncThreshold: The projection truncation threshold (default: 1.0e+1). See HuzinagaFDEProjectionPotential.h.
   *        - distantKinFunc: A flag whether not projected subsystems should be treated with a non additive kin. energy func. See HuzinagaFDEProjectionPotential.h.
   *        - useConvAcceleration: Turn the convergence acceleration (DIIS/Damping) on (default false). See FaTConvergenceAccerlerator.h.
   *        - diisStart: Density RMSD threshold for the start of the DIIS (default 5.0e-5).
   *        - diisEnd: Density RMSD threshold for the end of the DIIS (default 1.0e-4).
   *        - levelShiftParameter: The level-shift parameter (default 1.0e+6).
   *        - basisFunctionRatio: The basis function shell ratio for hybrid projection methods.
   *        - borderAtomThreshold: The threshold for environment orbitals on not distant atoms to consider them not distant.
   */
  FreezeAndThawTaskSettings settings;
private:
  std::vector<std::shared_ptr<SystemController> > _activeSystems;
  std::vector<std::shared_ptr<SystemController> > _passiveSystems;
};

} /* namespace Serenity */

#endif /* CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_ */
