/**
 * @file TopDownStaticEmbeddingTask.h
 *
 * @date Mar 21, 2024
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

#ifndef SERENITY_TOPDOWNSTATICEMBEDDINGTASK_H
#define SERENITY_TOPDOWNSTATICEMBEDDINGTASK_H

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/EmbeddingSettings.h"
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/SystemAdditionTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {

using namespace Serenity::Reflection;

struct TopDownStaticEmbeddingTaskSettings {
  TopDownStaticEmbeddingTaskSettings()
    : maxCycles(100), normThreshold(1e-5), writePairEnergies(false), useQuasiRestrictedOrbitals(true) {
    lcSettings.enforceHFFockian = true;
    loc.splitValenceAndCore = true;
    loc.useEnergyCutOff = false;
    lcSettings.linearScalingSigmaVector = true;
    lcSettings.method = Options::PNO_METHOD::NONE;
    lcSettings.embeddingSettings.embeddingMode = Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA;
    add.addOccupiedOrbitals = false;
  }
  REFLECTABLE((unsigned int)maxCycles, (double)normThreshold, (bool)writePairEnergies, (bool)useQuasiRestrictedOrbitals)
 public:
  LocalCorrelationSettings lcSettings;
  LocalizationTaskSettings loc;
  SystemSplittingTaskSettings split;
  SystemAdditionTaskSettings add;
};

class SystemController;

/**
 * @class TopDownStaticEmbeddingTask TopDownStaticEmbeddingTask.h
 * @tparam SCFMode Restricted vs. unrestricted.
 * @brief This task runs a supersystem DFT/HF calculation, selects orbital spaces and then performs an embedded
 *        calculation for the orbital spaces.
 */
template<Options::SCF_MODES SCFMode>
class TopDownStaticEmbeddingTask : public Task {
 public:
  /**
   * @brief Constructor.
   *
   * @param activeSystems          A list of all the active systems.
   * @param passiveSystems         A list of all the passive systems (never turning active).
   * @param supersystem            The (optional) supersystem.
   */
  TopDownStaticEmbeddingTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                             const std::vector<std::shared_ptr<SystemController>>& passiveSystems = {},
                             std::shared_ptr<SystemController> supersystem = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~TopDownStaticEmbeddingTask() = default;
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
  void visit(TopDownStaticEmbeddingTaskSettings& c, set_visitor v, std::string blockname) {
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
    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in TopDownStaticEmbeddingTaskSettings: " + blockname);
  }

  /**
   * @brief The settings.
   * - maxCycles           Maximum number of cycles for the resiudal equations.
   * - normThreshold       Convergence threshold for the residual equations (energy and abs. max. residual)
   * - writePairEnergies   (Coupled cluster only) Write the pair energies to file.
   * - useQuasiRestrictedOrbitals Use quasi restricted orbitals.
   * - lcSettings          The local correlation settings.
   * - loc                 The localization task settings.
   * - split               Settings to control the orbital partitioning.
   * - add                 Settings to control the supersystem construction.
   */
  TopDownStaticEmbeddingTaskSettings settings;

  /**
   * @brief Getter for the final energy.
   * @return The final energy.
   */
  double getFinalEnergy();

  /**
   * @brief Construct subsystems from a supersystem.
   * @param supersystem The supersystem.
   * @param activeSystems The active subsystems.
   * @param environmentSystems The environment subsystems.
   * @param loc The localization task settings.
   * @param split The system splitting task settings.
   * @param add The system addition task settings.
   * @param trunc The basis set truncation task settings.
   * @param useQuasiRestrictedOrbitals If true, quasi restricted orbitals are constructed.
   */
  static void setUpSubsystems(std::shared_ptr<SystemController> supersystem,
                              std::vector<std::shared_ptr<SystemController>> activeSystems,
                              std::vector<std::shared_ptr<SystemController>> environmentSystems,
                              const LocalizationTaskSettings& loc, const SystemSplittingTaskSettings& split,
                              const SystemAdditionTaskSettings& add, const BasisSetTruncationTaskSettings& trunc,
                              bool useQuasiRestrictedOrbitals = false);

 private:
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  std::vector<std::shared_ptr<SystemController>> _passiveSystems;
  std::vector<std::shared_ptr<SystemController>> _allSubsystems;
  const unsigned int _referenceIndex = 0;
  std::shared_ptr<SystemController> _supersystem;
  std::unique_ptr<double> _finalTotalEnergy;

  Eigen::VectorXd runCorrelationCalculations();
  void printResults(const Eigen::VectorXd& correlationEnergyContributions);
  void calculateReferenceEnergy();
  void updateSubsystemEnergy(unsigned int subsystemIndex, std::vector<std::shared_ptr<SystemController>> allSubsystems);
  void updateFrozenEnvironmentEnergies();
  static void cacheHCoreIntegrals(std::shared_ptr<SystemController> supersystem,
                                  std::vector<std::shared_ptr<SystemController>> activeSystems,
                                  std::vector<std::shared_ptr<SystemController>> environmentSystems);
};

} /* namespace Serenity */

#endif // SERENITY_TOPDOWNSTATICEMBEDDINGTASK_H
