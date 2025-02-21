/**
 * @file WavefunctionEmbeddingTask.h
 *
 * @author Moritz Bensberg
 * @date Jul 8, 2020
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

#ifndef TASKS_WAVEFUNCTIONEMBEDDINGTASK_H_
#define TASKS_WAVEFUNCTIONEMBEDDINGTASK_H_
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
#include "postHF/LocalCorrelation/LocalCorrelationController.h"
#include "settings/Reflection.h"
#include "tasks/LocalizationTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Serenity {

class SystemController;
class OrbitalPair;
class OrbitalTriple;

struct WavefunctionEmbeddingTaskSettings {
  WavefunctionEmbeddingTaskSettings()
    : normThreshold(1e-5),
      maxCycles(100),
      fullDecomposition(false),
      fromFragments(false),
      accurateInteraction(true),
      writePairEnergies(false) {
    loc.splitValenceAndCore = true;
  }
  REFLECTABLE((double)normThreshold, (unsigned int)maxCycles, (bool)fullDecomposition, (bool)fromFragments,
              (bool)accurateInteraction, (bool)writePairEnergies)
 public:
  std::vector<LocalCorrelationSettings> lcSettings;
  SystemSplittingTaskSettings split;
  LocalizationTaskSettings loc;

  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("WFEMB")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};
/**
 * @class
 * @brief A task that allows to perform DLPNO-in-DLPNO calculations using
 *        either DLPNO-CCSD or DLPNO-CCSD(T0). All thresholds may be adjusted
 *        for the specific subsystem and are then applied to the corresponding
 *        orbital singles, pairs and triples.
 */
class WavefunctionEmbeddingTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param supersystem The supersystem controller.
   * @param systems The subsystem controller.
   *
   * At the current moment both the supersystem and the subsytem controller have to be supplied.
   * However, one or the other may be dummies.
   */
  WavefunctionEmbeddingTask(std::shared_ptr<SystemController> supersystem,
                            std::vector<std::shared_ptr<SystemController>> systems);
  /**
   * @brief Default destructor.
   */
  ~WavefunctionEmbeddingTask();
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   *
   * Resolve setting list given by multiple blocks of the same type.
   */
  void visit(WavefunctionEmbeddingTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.split.visitAsBlockSettings(v, blockname))
      return;
    if (c.loc.visitAsBlockSettings(v, blockname))
      return;

    int lastCharIndex = blockname.length() - 1;
    int regionIndex = -1;
    if (lastCharIndex > 0)
      regionIndex = std::stoi(blockname.substr(lastCharIndex, 1));
    if (regionIndex > (int)c.lcSettings.size() || regionIndex < 0) {
      throw SerenityError((std::string) "Unknown block in WavefunctionEmbeddingTaskSettings: " + blockname);
    }
    std::string reducedBlockName = blockname.substr(0, 2);
    LocalCorrelationSettings& lcSettings = c.lcSettings[regionIndex];
    if (lcSettings.visitAsBlockSettings(v, reducedBlockName))
      return;

    // If reached, the keyword is unknown.
    throw SerenityError((std::string) "Unknown block in WavefunctionEmbeddingTaskSettings: " + blockname);
  }
  /**
   * @brief Settings.
   *   - normThreshold          Convergence threshold for the residual equations.
   *   - maxCycles              Maximum number of cycles for the residual equations.
   *   - fullDecomposition      Perform a full decomposition of the energy including HF/DFT contributions.
   *   - fromFragments          Ignore the supersystem and overwrite it with the union of the subsystems.
   *   - accurateInteraction    If true, the settings with the lower subsystem index are used for cross system pairs
   *                            otherwise the settings with the higher subsystem index are used.
   *   - lcSettings             The local correlation settings for each subsystem. It is assumed that a lower subsystem
   *                            index corresponds to more accurate settins.
   *   - split                  Settings for the system splitting task @see SystemSplittingTask.h.
   *   - loc                    Settings for the localization task @see LocalizationTask.h.
   */
  WavefunctionEmbeddingTaskSettings settings;
  /**
   * @brief Getter for the fragment-wise correlation energy.
   * @return The correlation energies.
   */
  const Eigen::VectorXd& getFragmentEnergies();
  /**
   * @brief Getter for the fragment-wise interaction-correlation energy.
   * @return The fragment-wise interaction-correlation energy.
   */
  const Eigen::VectorXd& getFragmentWiseInteractionEnergy();
  /**
   * @brief Getter for the total energy of the system (HF+correlation).
   *
   * @return The total energy.
   */
  double getTotalEnergy();

  std::shared_ptr<LocalCorrelationController> getLocalCorrelationController();

  //  const Eigen::SparseMatrix<int>& getSuperToSubsystemOccSortingMatrices();

 private:
  // The supersystem controller.
  std::shared_ptr<SystemController> _supersystem;
  // The subsystem controller.
  std::vector<std::shared_ptr<SystemController>> _systems;
  // The orbital to subsystem index map.
  Eigen::SparseMatrix<int> _orbitalIndexMap;
  // The orbital-wise subsystem assignments.
  Eigen::VectorXi _orbitalAssignments;
  /**
   * @brief Set up the subsystems.
   * @return The orbital to subsystem assignment.
   */
  Eigen::SparseMatrix<int> setUpSubsystems();
  /**
   * @brief Set up the supersystem.
   * @return The orbital to subsystem assignment.
   */
  Eigen::SparseMatrix<int> setUpSupersystem();
  /**
   * @brief Build the subsystem assignments for each orbital.
   * @return The subsystem assignments.
   */
  Eigen::VectorXi buildOrbitalAssignments();
  // The subsystem to orbital pair map.
  std::shared_ptr<Matrix<std::vector<std::shared_ptr<OrbitalPair>>>> _pairMatrix;
  // The list of all orbital pairs.
  std::vector<std::shared_ptr<OrbitalPair>> _allPairs;
  // The list of all orbital triples over multiple subsystems.
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<OrbitalTriple>>>> _crossTriples;
  // The list of all orbital triples over only one subsystem.
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<OrbitalTriple>>>> _diagonalTriples;
  // The auxiliary function selection thresholds for each atom.
  Eigen::VectorXd _orbitalWiseMullikenThresholds;
  // The shell selection thresholds for each orbital.
  Eigen::VectorXd _orbitalToShellThresholds;
  // The orbital-wise PAO thresholds.
  Eigen::VectorXd _orbitalWiseDOIPAOThresholds;
  // The energies of each fragment.
  Eigen::VectorXd _fragmentEnergies;
  // The interaction energies for each fragment.
  Eigen::VectorXd _interactionEnergies;
  /**
   * @brief Build the orbital pairs with their specific thresholds.
   * @return The orbital pair list.
   */
  std::vector<std::shared_ptr<OrbitalPair>> buildOrbitalPairs(unsigned int iSys, unsigned int jSys);
  /**
   * @brief Assign the orbital triples to the triples lists above.
   * @param localCorrelationController The central local correlation controller.
   */
  void sortOrbitalTriples(std::shared_ptr<LocalCorrelationController> localCorrelationController);
  /**
   * @brief Calculate the triples correction.
   * @param localCorrelationController  The central local correlation controller.
   * @return The triples correction.
   */
  double calculateTriplesCorrection(std::shared_ptr<LocalCorrelationController> localCorrelationController);
  /**
   * @brief Calculate the full energy decomposition for the given subsystem.
   * @param iSys The subsystem index.
   */
  void performFullDecomposition(unsigned int iSys);
  /**
   * @brief Construct the integral threshold vectors (see above) from the orbital-wise subsystem assignements.
   */
  void buildIntegralThresholdVectors();

  std::unique_ptr<double> _totalEnergy;

  std::shared_ptr<LocalCorrelationController> _localCorrelationController;

  void checkForMP2Run();
  bool _useMP2Calculator = false;
  bool _useCCCalculator = false;
  bool _onlyNone = false;
};

} /* namespace Serenity */

#endif /* TASKS_WAVEFUNCTIONEMBEDDINGTASK_H_ */
