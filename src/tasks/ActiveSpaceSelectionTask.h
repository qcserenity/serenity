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
      noIsolatedMOs(false),
      clusterThreshold(0.3),
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
              (Options::POPULATION_ANALYSIS_ALGORITHMS)populationAlgorithm, (bool)load, (bool)noIsolatedMOs,
              (double)clusterThreshold, (bool)checkDegeneracies, (double)degeneracyFactor, (bool)alignPiOrbitals,
              (unsigned int)alignCycles, (bool)usePiBias, (double)biasThreshold, (double)biasAverage,
              (unsigned int)alignExponent, (bool)kineticAlign, (bool)skipLocalization, (bool)splitValenceAndCore)
};

/**
 * @class ActiveSpaceSelectionTask ActiveSpaceSelectionTask.h
 * @brief Tries to select an active space from localized orbitals of structures along a reaction coordinate.
 *        The systems are written to scratch as <systemName>_Act and <systemName>_Env.
 *        The *_Env systems contain the orbitals which are found in all given systems. The *_Act systems the others.\n\n
 *
 *        MB: I tried various other criteria: DOI, naddXC, differences in XC-Fock matrix elements, LMP2 pair energy
 *            differences, dipole approximations to the LMP2 pair energy difference.\n
 *            Only the criteria which are still implemented worked reasonably well.\n\n
 *
 *        References:\n
 *            J. Chem. Phys. 2019, 150, 214106.
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
   * for actual sensible embedding calculations. populationAlgorithm             The algorithm used for the population
   * analysis. load                            Load systems from file. noIsolatedMOs                   Try to enforce
   * that no isolated environment orbitals appear. clusterThreshold                Population threshold that assigns
   * atoms to orbitals. Used for noIsolatedMOs. checkDegeneracies               Check for orbitals that are very similar
   * with respect to the comparison criteria. degeneracyFactor                Threshold scaling for degeneracy.
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

  /// @brief The matrices of the populations of the occupied orbitals on the atoms of the system.
  std::vector<SPMatrix<SCFMode>> _orbitalPopulations;
  /// @brief The kinetic energies of the occupied orbitals for the different systems.
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _kineticEnergies;
  /// @brief Final list of orbitals which do not have a partner in the other structures.
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> _unpairedOrbitals;
  /// @brief Mapping between the orbitals of the different systems.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> _completeOrbitalMap;
  /// @brief Final occupations for the active and environment systems. This has to be identical for every systems!
  /// (Sanity checks)
  std::shared_ptr<std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>> _finalOccupation;
  /// @brief The subsystem pairs. Only filled if "keepSystemPairs=true".
  std::vector<std::pair<std::shared_ptr<SystemController>, std::shared_ptr<SystemController>>> _systemPairs;

  /**
   * @brief Returns the access of alpha orbitals.
   * @param nOcc The occupations.
   * @return The access of alpha orbitals.
   */
  int getSpin(SpinPolarizedData<SCFMode, unsigned int> nOcc);
  /**
   * @brief Calculates the occupations of the subsystems from the orbital selection.
   * @param partnerVector The selection of active orbitals.
   * @return The occupation numbers for environment and active system (act : env)
   */
  std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
  getOrbitalOccupations(SpinPolarizedData<SCFMode, Eigen::VectorXi>& partnerVector);
  /**
   * @brief Calculates the orbital populations on the atoms
   *        of the systems and stores them in _orbitalPopulations for every system.
   */
  void buildOrbitalPopulations();
  /**
   * @brief Calculate, localize and align the orbitals.
   */
  void prepareOrbitals();
  /**
   * @brief Builds the map between orbitals of the subsystems from the selected criteria (localization/kinetic energy)
   */
  void buildOrbitalMap();
  /**
   *  @brief Builds a vector for every systems which shows which orbital has a partner in the other systems.
   */
  void reduceOrbitalMap();
  /**
   *  @brief Calculates the kinetic energy for every orbital of every systems and stores it in _kineticEnergies
   */
  void calculateKineticEnergy();
  /**
   * @brief Pre align orbitals and localize them.
   */
  void alignPiOrbitals();
  /**
   * @brief Generates the comparison criteria.
   */
  void produceComparisonCriteria();
  /**
   * @brief Does the orbital comparison step.
   */
  void compareOrbitals();
  /**
   * @brief Asserts that the orbtial occupations is identical
   *        nOccPair for every active and environment system over the course of the reaction.
   * @param nOccPair The occupation numbers which will be checked (act : env) against _finalOccupation
   */
  void assertOccupationConsistency(std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>> nOccPair);
  /**
   * @brief Prints information about the system selection.
   * @param system The system for which the information is printed.
   */
  void printSystemInformations(std::shared_ptr<SystemController> system);
  /**
   * @brief Ensures consistency of the orbital selection by checking the map
   *        of environment orbitals between the subsystems. Furthermore, the
   *        consistency of degeneracies between the orbitals are checked.
   */
  void backMapping();
  /**
   * @brief Checks if environment orbitals are isolated. If so, they are included
   *        in the active system. Isolation is determined via assigning atoms to the
   *        respective orbitals. If an environmental orbital has no atom in common with
   *        any other environment orbital, it is considered to be isolated.
   */
  void orbitalClustering();
  /**
   * @brief Initialize subsystems for a given supersystem.
   */
  std::shared_ptr<SystemController> initializeSubsystem(std::shared_ptr<SystemController> supersystemController,
                                                        std::shared_ptr<Geometry> geometry, bool active, int charge, int spin);
};

} /* namespace Serenity */

#endif /* TASKS_ACTIVESPACESELECTIONTASK_H_ */
