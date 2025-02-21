/**
 * @file GeneralizedDOSTask.h
 *
 * @author Moritz Bensberg
 * @date Sep 21, 2020
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

#ifndef TASKS_GENERALIZEDDOSTASK_H_
#define TASKS_GENERALIZEDDOSTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/SPMatrix.h"       //Orbital populations.
#include "math/Matrix.h"                  //Storage of fragments.
#include "settings/LocalizationOptions.h" //POPULATION_ANALYSIS_ALGORITHMS
#include "settings/Reflection.h"          //Reflections
#include "tasks/Task.h"                   //Base class.
/* Include Std and External Headers */
#include <Eigen/SparseCore>

namespace Serenity {
using namespace Serenity::Reflection;
/* Forward declarations */
class SystemController;
template<Options::SCF_MODES SCFMode>
class DirectOrbitalSelection;
class DOSOrbitalGroup;

struct GeneralizedDOSTaskSettings {
  GeneralizedDOSTaskSettings()
    : similarityLocThreshold({5e-2}),
      similarityKinEnergyThreshold({5e-2}),
      prioFirst(false),
      localizationThreshold(0.8),
      populationAlgorithm(Options::POPULATION_ANALYSIS_ALGORITHMS::IAOShell),
      checkDegeneracies(true),
      degeneracyFactor(1.0),
      usePiBias(true),
      biasThreshold(0.01),
      biasAverage(12.0),
      writeScores(false),
      scoreStart(1e-1),
      scoreEnd(1e-4),
      nTest(1000),
      mapVirtuals(false),
      writeGroupsToFile(false),
      bestMatchMapping(false) {
  }
  REFLECTABLE((std::vector<double>)similarityLocThreshold, (std::vector<double>)similarityKinEnergyThreshold,
              (bool)prioFirst, (double)localizationThreshold, (Options::POPULATION_ANALYSIS_ALGORITHMS)populationAlgorithm,
              (bool)checkDegeneracies, (double)degeneracyFactor, (bool)usePiBias, (double)biasThreshold,
              (double)biasAverage, (bool)writeScores, (double)scoreStart, (double)scoreEnd, (unsigned int)nTest,
              (bool)mapVirtuals, (bool)writeGroupsToFile, (bool)bestMatchMapping)
 public:
  /**
   * @brief Parse the settings from the input an instance of this class.
   * @param c The settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("DOS")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};
/**
 * @class
 * @brief A class that implements the DOS-procedure for any number of supersystems that are partitioned
 *        into at least two fragments each.
 *
 *        Since a supersystem may be partitioned into more than two
 *        fragments, this is the generalized DOS-procedure in contrast to the DOS procedure implemented in
 *        ActiveSpaceSelectionTask.\n\n
 *
 *        The generalized DOS (GDOS) runs DOS selections with increasingly tighter selection thresholds.
 *        The fragments are constructed from the additionally select orbitals in each iteration.
 */
template<Options::SCF_MODES SCFMode>
class GeneralizedDOSTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param acitveSystems         The supersystems.
   * @param environmentSystems    The ordered fragments. The first n fragments correspond to the first supersystem etc.,
   *                              where n is the number of fragments each supersystem is partitionied in.
   */
  GeneralizedDOSTask(std::vector<std::shared_ptr<SystemController>> activeSystems,
                     std::vector<std::shared_ptr<SystemController>> environmentSystems);
  /**
   * @brief Default destructor.
   */
  ~GeneralizedDOSTask();
  /**
   * @brief Run the task.
   */
  void run();
  /**
   * @brief Settings.
   * @param similarityLocThreshold       Localization threshold list.
   * @param similarityKinEnergyThreshold Kinetic energy threshold list.
   * @param prioFirst                    Prioritize the first fragment in the
   *                                     atom assignment after orbital selection.
   * @param localizationThreshold        If the occupied orbital population of the first
   *                                     fragment on a given atom exceeds this threshold
   *                                     the atom is assigned to the system (prioFirst=true).
   * @param populationAlgorithm          The algorithm used for the population analysis.
   * @param checkDegeneracies            Check for orbitals that are very similar
   *                                     with respect to the comparison criteria.
   * @param degeneracyFactor             Threshold scaling for the degeneracy check.
   * @param usePiBias                    Scale comparison threshold based on number of
   *                                     significant shells.
   * @param biasThreshold                SThreshold for the determination of an important shell.
   * @param biasAverage                  Scaling parameter for usePiBias.
   * @param writeScores                  If true, the scores at which each orbital is selected is written to file.
   *                                     The orbitals are not partitioned into subsystems. This is achieved with a
   *                                     large number of on a logarithmic scale tightly packed DOS thresholds.
   * @param scoreStart                   The start of the DOS-thresholds for the scan.
   * @param scoreEnd                     The end of the DOS-thresholds for the scan.
   * @param nTest                        The number of thresholds used in the scan.
   * @param mapVirtuals                  If true, the virtual orbitals are considered in the orbital mapping. By default
   *                                     false.
   * @param writeGroupsToFile            If true, a file is created containing the orbital set map between structures.
   * By default false.
   * @param bestMatchMapping             If true, the selection thresholds are optimized to provide a qualitative
   * orbital map, i.e., the thresholds are chosen such that they minimize the number of unmappable orbitals under the
   * constraint that they are smaller than "scoreStart". By default false.
   */
  GeneralizedDOSTaskSettings settings;

  const SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>>& getOrbitalGroups();
  const SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>>& getUnmappableOrbitalGroups();

  const std::vector<SpinPolarizedData<SCFMode, Eigen::SparseMatrix<int>>>& getSuperToSubsystemOccSortingMatrices();

 private:
  ///@brief The supersystems
  std::vector<std::shared_ptr<SystemController>> _supersystems;
  ///@brief The list of all fragments.
  std::vector<std::shared_ptr<SystemController>> _allFragments;
  ///@brief All Fragments ordered in a matrix. Row-index supersystem, col-index selection index.
  Matrix<std::shared_ptr<SystemController>> _fragments;
  ///@brief Number of fragments each supersystem is partitionined in.
  unsigned int _nFragments = 0;
  /**
   * @brief Construct the a _fragments-like matrix.
   * @param allFragments     The list of all fragments.
   * @param nSupersystems    The number of supersystems.
   * @param nFragmentsEach   The number of fragments for each supersystem.
   * @return The controller mapping matrix.
   */
  Matrix<std::shared_ptr<SystemController>>
  assignFragmentsToSupersystem(std::vector<std::shared_ptr<SystemController>> allFragments, unsigned int nSupersystems,
                               unsigned int nFragmentsEach);
  /**
   * @brief Check the given input and settings.
   */
  void checkInput();
  /**
   * @brief Calculate the orbital kinetic energies.
   * @param supersystems The supersystems.
   * @param nOrbitals    The number of orbitals to calculate the kinetic energy for.
   * @return The orbital kinetic energies.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>>
  calculateKineticEnergies(std::vector<std::shared_ptr<SystemController>> supersystems,
                           SpinPolarizedData<SCFMode, unsigned int> nOrbitals);
  /**
   * @brief Calculate the orbital-wise populations.
   * @param supersystems The supersystems.
   * @param alg          The population analysis algorithm used.
   * @return The orbital wise populations.
   */
  std::vector<SPMatrix<SCFMode>> calculateOrbitalPopulations(std::vector<std::shared_ptr<SystemController>> supersystems,
                                                             Options::POPULATION_ANALYSIS_ALGORITHMS alg);
  /**
   * @brief Assign the index of the current iteration to the orbitals which were selected in this iteration.
   * @param finalAssignment The final assignments to manipulate.
   * @param newAssignment   The current DOS selection.
   * @param assignIndex     The index of the current DOS iteration.
   */
  void setAssignIndices(std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& finalAssignment,
                        const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& newAssignment, int assignIndex);
  /**
   * @brief Split the supersystem geometries based on the assignment.
   * @param assignments The assignemtns.
   */
  void splitGeometryIntoFragments(const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments);
  /**
   * @brief Split the occupied orbitals based on the assignemt.
   * @param assignments The assignments.
   */
  void splitElectronicStructureIntoFragments(const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments);
  /**
   * @brief Check if the occupations of all fragments is consistent.
   * @return True if no errors are present, else false.
   */
  bool checkFinalOccupations();
  /**
   * @brief Run the iterative selection.
   * @return The final assignments.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
  getIterativeAssignments(std::vector<double> locs, std::vector<double> kins, DirectOrbitalSelection<SCFMode>& dos);
  /**
   * @brief Getter for the number of orbitals for which populations were calculated.
   * @param populations The populations.
   * @return The number of orbitals for which populations were calculated.
   */
  SpinPolarizedData<SCFMode, unsigned int> getMaxNumberOfOrbitalsFromPopulations(const std::vector<SPMatrix<SCFMode>>& populations);
  /*
   * Numerical scan of the GDOS-function. The result is written to file.
   */
  void writeScoresToFile(DirectOrbitalSelection<SCFMode>& dos);
  void writeFile(SPMatrix<SCFMode> scores);
  void writeMatrix(Eigen::MatrixXd scores, std::string fileName = "");

  void writeGroupsToFile();
  void writeGroupSetToFile(std::string fileName, std::vector<std::shared_ptr<DOSOrbitalGroup>> groups);

  SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>> _orbitalGroups;
  SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>> _unmappableOrbitalGroups;

  std::unique_ptr<std::vector<SpinPolarizedData<SCFMode, Eigen::SparseMatrix<int>>>> _superToSubsystemOccSortingMatrices;

  std::unique_ptr<std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>> _finalAssignments;
};

} /* namespace Serenity */

#endif /* TASKS_GENERALIZEDDOSTASK_H_ */
