/**
 * @file ActiveSpaceSelectionTask.h
 *
 * @date Sep 11, 2018
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

#ifndef TASKS_ACTIVESPACESELECTIONTASK_H_
#define TASKS_ACTIVESPACESELECTIONTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/SPMatrix.h"
#include "settings/LocalizationOptions.h"
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
template<Options::SCF_MODES T>
class ElectronicStructure;

struct ActiveSpaceSelectionTaskSettings {
  ActiveSpaceSelectionTaskSettings()
    : similarityLocThreshold(5e-2),
      similarityKinEnergyThreshold(5e-2),
      locType(Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO),
      localizationThreshold(0.8),
      populationAlgorithm(Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell),
      load(false),
      checkDegeneracies(false),
      degeneracyFactor(10.0),
      alignPiOrbitals(false),
      alignCycles(10),
      usePiBias(false),
      biasThreshold(0.01),
      biasAverage(12.0),
      alignExponent(2),
      kineticAlign(false),
      skipLocalization(false),
      splitValenceAndCore(false) {
  }
  REFLECTABLE((double)similarityLocThreshold, (double)similarityKinEnergyThreshold,
              (Options::ORBITAL_LOCALIZATION_ALGORITHMS)locType, (double)localizationThreshold,
              (Options::POPULATION_ANALYSIS_ALGORITHMS)populationAlgorithm, (bool)load, (bool)checkDegeneracies,
              (double)degeneracyFactor, (bool)alignPiOrbitals, (unsigned int)alignCycles, (bool)usePiBias,
              (double)biasThreshold, (double)biasAverage, (unsigned int)alignExponent, (bool)kineticAlign,
              (bool)skipLocalization, (bool)splitValenceAndCore)
};

/**
 * @class ActiveSpaceSelectionTask ActiveSpaceSelectionTask.h
 * @brief Tries to select an active space from localized orbitals of structures along a reaction coordinate.
 *
 *        The systems are written to scratch as <systemName>_Act and <systemName>_Env if no subsystems are provided..
 *        The *_Env systems contain the orbitals which are found in all given systems. The *_Act systems the others.\n\n
 *
 *        See DirectOrbitalSelection and GeneralizedDOSTask for more information.
 */
template<Options::SCF_MODES SCFMode>
class ActiveSpaceSelectionTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param supersystems        The systems which will be analyzed in order to select an active space for them.
   *                            At least two systems are needed!
   * @param activeSystems       The system controllers to which the active systems will be assigned.
   *                            If none given, a system controller will be constructed and its electronic-structure
   *                            printed to disk. Names will be supersystem-name + "_act".
   * @param environmentSystems  The system controllers to which the environment systems will be assigned.
   *                            If none given, a system controller will be constructed and its electronic-structure
   *                            printed to disk. Names will be supersystem-name + "_env".
   */
  ActiveSpaceSelectionTask(std::vector<std::shared_ptr<SystemController>> supersystems,
                           std::vector<std::shared_ptr<SystemController>> acitveSystems,
                           std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Default Destructor.
   */
  virtual ~ActiveSpaceSelectionTask() = default;
  /**
   * @brief Execute the task.
   */
  void run();
  /**
   * @brief Settings.
   *  similarityLocThreshold          Threshold for the difference in orbital localization between occupied orbitals.
   *  similarityKinEnergyThreshold    Threshold for the difference in kinetic energy.
   *  locType                         Localization type. See Options.h for options.
   *  localizationThreshold           Threshold for the assignment of atoms to the subsystems. This is purely cosmetic
   * for actual sensible embedding calculations.
   *  populationAlgorithm             The algorithm used for the population
   * analysis.
   *  load                            Load systems from file.
   *  checkDegeneracies               Check for orbitals that are very similar
   * with respect to the comparison criteria.
   *  degeneracyFactor                Threshold scaling for degeneracy.
   *  alignPiOrbitals                 Pre-align orbitals.
   *  alignCycles                     Number of align iterations before the procedure cancels if it does not converge.
   *  usePiBias                       Scale comparison threshold based on number of significant shells.
   *  biasThreshold                   Threshold for the determination of an important shell.
   *  biasAverage                     Scaling parameter for usePiBias.
   *  alignExponent                   Exponent used in orbital alignment. See LocalizationTaskSettings.
   *  skipLocalization                The orbitals of the systems are used with SCF, alignment or localization.
   *  splitValenceAndCore             Split valence and core orbitals during alignment AND localization.
   */
  ActiveSpaceSelectionTaskSettings settings;

  /// @brief A flag to keep the subsystem pairs.
  bool keepSystemPairs = false;
  /**
   * @brief Getter for the system pairs.
   * @return The system pairs.
   */
  std::vector<std::pair<std::shared_ptr<SystemController>, std::shared_ptr<SystemController>>> getSystemPairs() {
    if (!keepSystemPairs)
      throw SerenityError("System pairs were not saved! Nothing to return here! Adjust the settings of the task!");
    return _systemPairs;
  }

 private:
  /// @brief The system controllers of the systems which are compared.
  std::vector<std::shared_ptr<SystemController>> _supersystems;
  /// @brief The system controllers to store the selection results in.
  std::vector<std::shared_ptr<SystemController>> _activeSystems;
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  /// @brief The subsystem pairs. Only filled if "keepSystemPairs=true".
  std::vector<std::pair<std::shared_ptr<SystemController>, std::shared_ptr<SystemController>>> _systemPairs;
  /**
   * @brief Calculate, localize and align the orbitals.
   */
  void prepareOrbitals();
  /**
   * @brief Pre align orbitals and localize them.
   */
  void alignPiOrbitals();
  /**
   * @brief Initialize subsystems for a given supersystem.
   */
  std::shared_ptr<SystemController> initializeSubsystem(std::shared_ptr<SystemController> supersystemController,
                                                        std::string namePostfix);
  /**
   * @brief Construct a SystemController-vector for the GDOS.
   *        If necessary, new subsystem controller are constructed.
   * @return The systemController-vector.
   */
  std::vector<std::shared_ptr<SystemController>> sortOrCreateAllSubsystemController();
};

} /* namespace Serenity */

#endif /* TASKS_ACTIVESPACESELECTIONTASK_H_ */
