/**
 * @file DFTEmbeddedLocalCorrelationTask.h
 *
 * @date Apr 4, 2021
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

#ifndef SERENITY_DFTEMBEDDEDLOCALCORRELATIONTASK_H
#define SERENITY_DFTEMBEDDEDLOCALCORRELATIONTASK_H

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/CorrelatedMethodsOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Reflection.h"
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/SystemAdditionTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;
struct DFTEmbeddedLocalCorrelationTaskSettings {
  DFTEmbeddedLocalCorrelationTaskSettings()
    : maxCycles(100), normThreshold(1e-5), writePairEnergies(false), runFaT(false), fromSupersystem(true) {
    lcSettings.enforceHFFockian = true;
    loc.splitValenceAndCore = true;
    lcSettings.linearScalingSigmaVector = true;
    lcSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
    add.addOccupiedOrbitals = false;
    trunc.truncAlgorithm = Options::BASIS_SET_TRUNCATION_ALGORITHMS::NONE;
  }
  REFLECTABLE((unsigned int)maxCycles, (double)normThreshold, (bool)writePairEnergies, (bool)runFaT, (bool)fromSupersystem)
 public:
  LocalCorrelationSettings lcSettings;
  LocalizationTaskSettings loc;
  SystemSplittingTaskSettings split;
  SystemAdditionTaskSettings add;
  BasisSetTruncationTaskSettings trunc;
};
/**
 * @class
 * @brief Perform a DFT-embedded local-correlation calculation. This task only servers as an input simplification
 *        by calling the appropriate other tasks.
 */
class DFTEmbeddedLocalCorrelationTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem        The system controller.
   * @param environmentSystems  The environment systems.
   * @param supersystem         The (optional) supersystem.
   * Assigns private variables only.
   */
  DFTEmbeddedLocalCorrelationTask(std::shared_ptr<SystemController> activeSystem,
                                  std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                  std::shared_ptr<SystemController> supersystem = nullptr);
  /**
   * @brief Default destructor
   */
  ~DFTEmbeddedLocalCorrelationTask();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(DFTEmbeddedLocalCorrelationTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.lcSettings.visitAsBlockSettings(v, blockname))
      return;
    if (c.loc.visitAsBlockSettings(v, blockname))
      return;
    if (c.add.visitAsBlockSettings(v, blockname))
      return;
    if (c.split.visitAsBlockSettings(v, blockname))
      return;
    if (c.trunc.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in DFTEmbeddedLocalCorrelationTaskSettings: " + blockname);
  }
  /**
   * @brief The task settings.
   * - maxCycles           Maximum number of cycles for the resiudal equations.
   * - normThreshold       Convergence threshold for the residual equations (energy and abs. max. residual)
   * - writePairEnergies   (Coupled cluster only) Write the pair energies to file.
   * - lcSettings          The local correlation settings.
   * - loc                 The localization task settings.
   * - split               The system splitting task settings.
   * - add                 The system addition task settings.
   * - trunc               The basis set truncation task settings.
   */
  DFTEmbeddedLocalCorrelationTaskSettings settings;
  /**
   * @see Task
   */
  void run();

 private:
  // The system controller.
  std::shared_ptr<SystemController> _activeSystem;
  // The environment system controllers.
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  // The supersystem controller. This will only be initialized if the subsystems are to be constructed
  // from a supersystem guess.
  std::shared_ptr<SystemController> _supersystem;
  // Set up subsystem molecular orbitals from a supersystem calculation.
  void setUpSubsystems();
};

} /* namespace Serenity */
#endif // SERENITY_DFTEMBEDDEDLOCALCORRELATIONTASK_H
