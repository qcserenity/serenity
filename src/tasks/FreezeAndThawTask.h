/**
 * @file FreezeAndThawTask.h
 *
 * @date Oct 15, 2015
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

#ifndef CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_
#define CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_

/* Include Serenity Internal Headers */
#include "dft/Functional.h"
#include "misc/SerenityError.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory> //smart ptr.
#include <vector>

namespace Serenity {
/* Forward declaration */
class SystemController;

using namespace Serenity::Reflection;

struct FreezeAndThawTaskSettings {
  FreezeAndThawTaskSettings()
    : maxCycles(50),
      convThresh(1.0e-6),
      gridCutOff(-1.0),
      smallSupersystemGrid(false),
      basisExtThresh(5.0e-2),
      extendBasis(false),
      useConvAcceleration(false),
      diisStart(5.0e-5),
      diisEnd(1e-4),
      calculateSolvationEnergy(false),
      calculateUnrelaxedMP2Density({}),
      mp2Type(Options::MP2_TYPES::LOCAL),
      keepCoulombCache(false),
      finalEnergyEvaluation(true),
      printResults(true),
      onlyFinalEnergyEvaluation(false) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::PW91;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NADD_FUNC;
    lcSettings.pnoSettings = Options::PNO_SETTINGS::TIGHT;
    lcSettings.method = Options::PNO_METHOD::DLPNO_MP2;
  }
  REFLECTABLE((unsigned int)maxCycles, (double)convThresh, (double)gridCutOff, (bool)smallSupersystemGrid,
              (double)basisExtThresh, (bool)extendBasis, (bool)useConvAcceleration, (double)diisStart, (double)diisEnd,
              (bool)calculateSolvationEnergy, (std::vector<bool>)calculateUnrelaxedMP2Density, (Options::MP2_TYPES)mp2Type,
              (bool)keepCoulombCache, (bool)finalEnergyEvaluation, (bool)printResults, (bool)onlyFinalEnergyEvaluation)
 public:
  EmbeddingSettings embedding;
  LocalCorrelationSettings lcSettings;
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
  FreezeAndThawTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                    const std::vector<std::shared_ptr<SystemController>>& passiveSystems = {});
  /**
   * @brief Default destructor.
   */
  virtual ~FreezeAndThawTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(FreezeAndThawTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    // If reached, the blockname is unknown.
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    throw SerenityError((std::string) "Unknown block in FreezeAndThawTaskSettings: " + blockname);
  }

  /**
   * @brief The settings/keywords for the FreezeAndThawTask:
   *        - maxCycles: The maximum number of FAT iterations (default: 50)
   *        - gridCutOff :  Modifies the super-system grid used to evaluate the FDE potential. See FDETask.
   *        - basisExtensionThreshold: The threshold for the basis set extension (default: 5.0 e-2). See BasisExtension.h
   *        - extendBasis: A flag whether the basis should be extended (default: false). See BasisExtension.h
   *        - truncateProjector: A flag whether the projector should be truncated (default: false). See
   * HuzinagaProjectionPotential.h.
   *        - projectionTruncThreshold: The projection truncation threshold (default: 1.0e+1). See
   * HuzinagaProjectionPotential.h.
   *        - distantKinFunc: A flag whether not projected subsystems should be treated with a non additive kin. energy
   * func. See HuzinagaProjectionPotential.h.
   *        - useConvAcceleration: Turn the convergence acceleration (DIIS/Damping) on (default false). See
   * FaTConvergenceAccerlerator.h.
   *        - diisStart: Density RMSD threshold for the start of the DIIS (default 5.0e-5).
   *        - diisEnd: Density RMSD threshold for the end of the DIIS (default 1.0e-4).
   *        - calculateSolvationEnergy: Calculate only the interaction of the first active system with the environment
   *                                    and the energy of the first active system.
   *        - keepCoulombCache: The Fock matrix contributions of the passive systems via their Coulomb interaction
   *                            is not deleted.
   *        - embedding: The embedding settings. See settings/EmbeddingSettings.h for details.
   *        - onlyFinalEnergyEvaluation : Run only the final energy evaluation.
   *        - printResults : Print the final results to the std-output.
   */
  FreezeAndThawTaskSettings settings;

 private:
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  std::vector<std::shared_ptr<SystemController>> _passiveSystems;
  /**
   * @brief Calculate the non-additive dispersion correction.
   *
   * This step can become computationally expensive for thousands of atoms.
   * Thus, we only want to do it once at the end of the freeze-and-thaw iterations.
   * Furthermore, this is different if the settings <calculateSolvationEnergy>
   * is chosen as true.
   */
  void calculateNonAdditiveDispersionCorrection();
  /**
   * @brief Clean up the passive Coulomb cache after finishing the iterations.
   */
  void cleanUp();
  /**
   * @brief Calculate the final energy.
   */
  void finalEnergyEvaluation();
};

} /* namespace Serenity */

#endif /* CONFIGURATION_TASKS_FREEZEANDTHAWTASK_H_ */
