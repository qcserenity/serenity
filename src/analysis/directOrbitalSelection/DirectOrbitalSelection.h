/**
 * @file DirectOrbitalSelection.h
 *
 * @author Moritz Bensberg
 * @date Sep 18, 2020
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

#ifndef DIRECTORBITALSELECTION_DIRECTORBITALSELECTION_H_
#define DIRECTORBITALSELECTION_DIRECTORBITALSELECTION_H_

/* Include Serenity Internal Headers */
#include "data/matrices/SPMatrix.h" //Orbital populations.
/* Include Std and External Headers */
#include <Eigen/SparseCore> // Sparse matrices.
#include <iomanip>          // setw
#include <vector>           // std::vector.

namespace Serenity {

/**
 * @class
 * @brief A group of orbital indices for multiple systems along a reaction coordinate. The orbitals corresponding to
 * these indices are considered indistinguishable by the DOS.
 */
class DOSOrbitalGroup {
 public:
  /**
   * @brief Constructor.
   * @param systemWiseOrbitalIndices The system-wise indices of the orbitals. The outer index must correspond to the
   * system. The inner index enumerates the orbitals. Each orbital set must have the same size.
   */
  DOSOrbitalGroup(std::vector<std::vector<unsigned int>> systemWiseOrbitalIndices);
  DOSOrbitalGroup() = default;

  /**
   * @brief Construct orbital groups from a vector that denotes the assignment of orbitals to a subsystem.
   * @param assignment  The orbital to subsystem assignment for each system.
   * @param targetIndex The subsystem target index, i.e., orbitals assigned to this subsystem are collected
   *                    into the orbital group for each system.
   * @return The orbital groups (all orbitals in the given subsystem).
   */
  template<Options::SCF_MODES SCFMode>
  static SpinPolarizedData<SCFMode, DOSOrbitalGroup>
  orbitalGroupFromAssignment(const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignment,
                             unsigned int targetIndex);
  /**
   * @brief Construct orbital groups from an orbital-to-orbital map between subsystems, the system-wise assignment and
   * the subsystem index to consider. This provides a more detailed splitting of the orbital groups than
   * orbitalGroupFromAssignment(...).
   * @param map         The orbital-to-orbital map between systems.
   * @param assignment  The orbital to subsystem assignment.
   * @param targetIndex The subsystem index.
   * @param The orbital groups.
   */
  template<Options::SCF_MODES SCFMode>
  static SpinPolarizedData<SCFMode, std::vector<std::shared_ptr<DOSOrbitalGroup>>>
  orbitalGroupsFromMapWithAssignment(const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
                                     const SpinPolarizedData<SCFMode, Eigen::VectorXi>& assignment, unsigned int targetIndex);
  /**
   * @brief Construct sparse sorting matrices between orbital indices and group indies.
   * @param orbitalGroups The orbital group list.
   * @param nOcc          The total number of orbitals.
   */
  static std::vector<Eigen::SparseMatrix<int>>
  buildOrbitalToOrbitalGroupSortingMatrix(std::vector<std::shared_ptr<DOSOrbitalGroup>> orbitalGroups, const unsigned int nOcc);
  /**
   * @brief Collect pair energies group wise.
   * @param pairEnergies All pair energies.
   * @return The sum of all pair energies for which at least one orbital index is in the given group.
   */
  Eigen::VectorXd getGroupEnergiesFromPairEnergies(const std::vector<Eigen::MatrixXd>& pairEnergies) const;

  /**
   * @brief Collect the change in pair energies for the given orbital group between all systems.
   * @param pairEnergies All pair energies.
   * @return A matrix that provides the change in group pair energies between systems.
   */
  Eigen::MatrixXd getGroupEnergyDifferenceMatrix(const std::vector<Eigen::MatrixXd>& pairEnergies) const;
  /**
   * @brief Getter for the orbital indices for the given system index in the group.
   * @param systemIndex The system index.
   * @return The orbital indices.
   */
  const std::vector<unsigned int>& getOrbitalIndicesForSystem(unsigned int systemIndex);
  /**
   * @brief Getter for the orbital indices for each system.
   * @return The system orbital indices. The outer index corresponds to the system.
   */
  const std::vector<std::vector<unsigned int>>& getSystemWiseIndices() const;
  /**
   * @brief Add a new set of orbitla indices for a new system.
   * @param systemIndices The orbital indices. The number of orbital indices must be the same as the
   *                      number of indices for all other sets in this group.
   */
  void addSystemWiseIndices(std::vector<unsigned int> systemIndices);
  /**
   * @brief Update the orbital indices according to the permutation maps provided.
   * @param newXorigMap The permutation maps.
   */
  void updateOrbitalIndices(std::vector<Eigen::SparseMatrix<int>> newXorigMap);
  /**
   * @brief Getter for the number of orbitals in this group.
   * @return The number of orbitals mapped through this group.
   */
  unsigned int getNOrbitals() const;
  /**
   * @brief Getter for the number of systems contributing to this group.
   * @return The number of systems.
   */
  unsigned int getNSystems() const;
  /**
   * @brief Split a group into virtual and occupied subgroups.
   * @param nOcc The index of the first virtual orbital / the number of occupied orbitals.
   * @return The virtual subgroup.
   */
  DOSOrbitalGroup splitOffVirtualOrbitals(const unsigned int nOcc);
  /**
   * @brief Output stream to print the orbital indices in this group.
   *
   * Format: orb-of-sys1 orb-of-sys2 orb-of-sys3 ...
   *         orb-of-sys1 orb-of-sys2 orb-of-sys3 ...
   *         orb-of-sys1 orb-of-sys2 orb-of-sys3 ...
   */
  friend std::ostream& operator<<(std::ostream& os, const DOSOrbitalGroup& group) {
    const auto& systemWiseOrbitalIndices = group.getSystemWiseIndices();
    for (unsigned int iOrb = 0; iOrb < group.getNOrbitals(); ++iOrb) {
      for (unsigned int iSys = 0; iSys < group.getNSystems(); ++iSys) {
        os << " " << systemWiseOrbitalIndices[iSys][iOrb] << std::setw(6);
      }
      os << "\n";
    }
    return os;
  }
  /**
   * @brief Collect the orbital indices for multiple groups. The first index corresponds to the group, the second to the
   * system, and the third enumerates the orbitals in the group.
   * @param The orbital groups.
   * @return All orbital indices in this group.
   */
  static std::vector<std::vector<std::vector<unsigned int>>>
  getIndicesFromGroups(std::vector<std::shared_ptr<DOSOrbitalGroup>> groups);

 private:
  std::vector<std::vector<unsigned int>> _systemWiseOrbitalIndices = {};
};

/**
 * @class
 * @brief Implements the DOS algorithm.
 *
 *        References:\n
 *            J. Chem. Phys. 2019, 150, 214106.\n
 *            J. Chem. Theory Comput. 2020, 16, 3607--3619\n\n
 *
 *        MB: I tried various other criteria: DOI, naddXC, differences in XC-Fock matrix elements, LMP2 pair energy
 *            differences, dipole approximations to the LMP2 pair energy difference.\n
 *            Only the criteria which are still implemented worked reasonably well.\n\n
 */
template<Options::SCF_MODES SCFMode>
class DirectOrbitalSelection {
 public:
  /**
   * @brief Constructor.
   * @param orbitalPopulations     The orbital-wise populations.
   * @param kineticEnergies        The orbital-wise kinetic energies.
   * @param nOccupiedOrbitals      The number of occupied orbitals/orbitals to compare.
   * @param usePiBias              Flag to enable threshold scaling by the number of important shells..
   * @param biasThreshold          Threshold to consider a shell as important for the orbital.
   * @param biasAverage            Threshold scaling parameter.
   * @param checkDegeneracies      Check the occupied orbitals for degeneracies.
   * @param degeneracyFactor       Scaling factor for the degeneracy check.
   * @param excludeOccVirtMappings Orbital mappings between virtual and occupied orbitals are not allowed.
   */
  DirectOrbitalSelection(const std::vector<SPMatrix<SCFMode>>& orbitalPopulations,
                         const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> kineticEnergies,
                         const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals, bool usePiBias,
                         double biasThreshold, double biasAverage, bool checkDegeneracies, double degeneracyFactor,
                         bool excludeOccVirtMappings = true);
  /**
   * @brief Default destructor.
   */
  ~DirectOrbitalSelection();
  /**
   * @brief Getter for the orbital assignment.
   * @param locThreshold  The orbital population threshold.
   * @param kinThreshold  The orbital-wise kinetic energy threshold.
   * @return The DOS-assignments. 1 -> not selected, 0 -> selected.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> getOrbitalAssignments(double locThreshold, double kinThreshold);
  /**
   * @brief Construct the assignment with the smallest number of unmappable orbitals. Note that this uses a binary
   * search and both selection thresholds are assumed to be equal.
   * @param threshold     The orbital-wise kinetic and localization threshold. Will be overwritten by the function.
   * @param map Surjective maps between the system orbitals.
   *            First two indices correspond to the subsystems, the last two
   *            (Eigen::MatrixXi-indices) to the orbitals.
   * @param minThreshold  Lower limit for the threshold.
   * @param maxThreshold  Upper limit for the threshold.
   * @return The best matching assignment.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
  getBestMatchAssignments(double& threshold, double& virtThreshold,
                          std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
                          const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals, double minThreshold = 5e-2,
                          double maxThreshold = 5e-1);
  /**
   * @brief Getter for the orbital assignment.
   * @param locThreshold  The orbital population threshold.
   * @param kinThreshold  The orbital-wise kinetic energy threshold.
   * @param map Surjective maps between the system orbitals.
   *            First two indices correspond to the subsystems, the last two
   *            (Eigen::MatrixXi-indices) to the orbitals.
   * @return The DOS-assignments. 1 -> not selected, 0 -> selected.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
  getOrbitalAssignments(double locThreshold, double kinThreshold,
                        std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map);

 private:
  /**
   * @brief Calculate pair-wise kinetic energy differences.
   * @param nOccupiedOrbitals The number of occuped orbitals.
   */
  void calculateKinScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals);
  /**
   * @brief Calculate pair-wise orbital population differences.
   * @param nOccupiedOrbitals The number of occuped orbitals.
   */
  void calculateLocScore(const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals);
  /**
   * @brief Builds the map between orbitals of the subsystems from the selected criteria (localization/kinetic energy)
   */
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> buildOrbitalMap(double locThreshold, double kinThreshold);
  /**
   *  @brief Builds a vector for every systems which shows which orbital has a partner in the other systems.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
  reduceOrbitalMap(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map);
  /**
   * @brief Ensures consistency of the orbital selection by checking the map
   *        of environment orbitals between the subsystems. Furthermore, the
   *        consistency of degeneracies between the orbitals are checked.
   */
  void backMapping(std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments,
                   std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> map);
  /**
   * @brief Construct the orbital sets which we intend to map from the orbital self-map (comparison
   * of the orbitals within the same system). Orbitals that are considered unmappable in the assignments
   * will never be included in such a set.
   * @param map          All orbital maps.
   * @param assignments  The orbital assignment.
   * @return The orbital sets as list of indices for each system. List ordering: System-index, set index, orbital index
   * in set.
   */
  std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>
  buildSystemWiseOrbitalSets(const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& map,
                             const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments);
  /**
   * @brief Construct the map between orbital sets.
   * @param orbitalMap   The original map between orbitals.
   * @param orbitalSets  The definition of orbital sets.
   * @return The maps between orbital sets.
   */
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>
  buildOrbitalSetMap(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalMap,
                     const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets);
  /**
   * @brief Check if there is a map between the sets iSet and jSet.
   * @param iSet The list of indices of the first set.
   * @param jSet The list of indices of the second set.
   * @param map  The map between the orbitals.
   * @return True, if all orbitals in the sets are mapped to all orbitals in the respective other set through the map.
   * Otherwise, false.
   */
  bool mapBetweenSets(const std::vector<unsigned int>& iSet, const std::vector<unsigned int>& jSet, const Eigen::MatrixXi& map);
  /**
   * @brief Ensure that the mapping between orbital sets includes only sets that can be mapped for all system
   * combinations.
   * @param orbitalSetMap The map between orbital sets.
   * @return The orbital set assignment, i.e, a list of indices that are 0 if the corresponding set cannot always be
   * mapped.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>
  enforceConsistentOrbitalSetMap(const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalSetMap);
  /**
   * @brief Translate the orbital set map back into an orbital map, i.e., create a map in which only maps between
   * orbitals are defined if there is also a map between the corresponding orbital sets and this set can be mapped
   * between all structures.
   * @param orbitalSetMap         The orbital set map.
   * @param orbitalSetAssignment  The orbital set assignment.
   * @param orbitalSets           The definition of the orbital sets.
   * @return The orbital map.
   */
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>
  setMapsToOrbitalMap(const std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>>& orbitalSetMap,
                      const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& orbitalSetAssignment,
                      const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets);
  /**
   * @brief Update the orbital assignments with the set assignments, i.e., if an orbital can be mapped through
   * the original orbital map but its orbital set cannot always be mapped, the orbital's assignment is changed to 0.
   * @param orbitalSets           The orbital sets.
   * @param orbitalSetAssignment  The orbital set assignment.
   * @param assignments           The original orbital assignment. Changed in place.
   */
  void updateOrbitalAssignmentWithSetAssignment(
      const std::vector<SpinPolarizedData<SCFMode, std::vector<std::vector<unsigned int>>>>& orbitalSets,
      const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& orbitalSetAssignment,
      std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>>& assignments);

  ///@brief Initialize the assignment vectors.
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> initializeAssignments();
  ///@brief Initialize the orbital map.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> initializeOrbitalMap();
  ///@brief Initialize a score field.
  void initializeScores(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>>& scores,
                        const SpinPolarizedData<SCFMode, unsigned int>& nOrbitals);
  SpinPolarizedData<SCFMode, unsigned int> initializeNOrbitals();
  ///@brief The orbital-wise populations.
  std::vector<SPMatrix<SCFMode>> _orbitalPopulations;
  ///@brief The orbital-wise kinetic energies.
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _kineticEnergies;
  ///@brief The population-difference score.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>> _locScore;
  ///@brief The kinetic-energy-difference score.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>> _kinScore;
  /// The DOS settings.
  const bool _usePiBias;
  const double _biasThreshold;
  const double _biasAverage;
  const bool _checkDegeneracies;
  const double _degeneracyFactor;
  const bool _excludeOccVirtMappings;
  const SpinPolarizedData<SCFMode, unsigned int> _nOrbitals;
};

} /* namespace Serenity */

#endif /* DIRECTORBITALSELECTION_DIRECTORBITALSELECTION_H_ */
