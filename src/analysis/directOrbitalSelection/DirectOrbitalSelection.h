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

namespace Serenity {
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
   * @param nOccupiedOrbitals      The number of occupied orbitals.
   * @param usePiBias              Flag to enable threshold scaling by the number of important shells..
   * @param biasThreshold          Threshold to consider shella as important for the orbital.
   * @param biasAverage            Threshold scaling parameter.
   * @param checkDegeneracies      Check the occupied orbitals for degeneracies.
   * @param degeneracyFactor       Scaling factor for the degeneracy check.
   */
  DirectOrbitalSelection(const std::vector<SPMatrix<SCFMode>>& orbitalPopulations,
                         const std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXd>> kineticEnergies,
                         const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals, bool usePiBias,
                         double biasThreshold, double biasAverage, bool checkDegeneracies, double degeneracyFactor);
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
   * @brief Getter for the orbital assignment.
   * @param locThreshold  The orbital population threshold.
   * @param kinThreshold  The orbital-wise kinetic energy threshold.
   * @param map Surjective maps between the system occupied orbitals.
   *            First two indices correspond to the subsystems, the last two
   *            (Eigen::MatrixXi-indices) to the occupied orbitals.
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
  ///@brief Initialize the assignment vectors.
  std::vector<SpinPolarizedData<SCFMode, Eigen::VectorXi>> initializeAssignments();
  ///@brief Initialize the orbital map.
  std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXi>>> initializeOrbitalMap();
  ///@brief Initialize a score field.
  void initializeScores(std::vector<std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>>& scores,
                        const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals);
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
};

} /* namespace Serenity */

#endif /* DIRECTORBITALSELECTION_DIRECTORBITALSELECTION_H_ */
