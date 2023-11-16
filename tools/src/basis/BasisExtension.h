/**
 * @file BasisExtension.h
 *
 * @date Jan 11, 2018
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

#ifndef BASIS_BASISEXTENSION_H_
#define BASIS_BASISEXTENSION_H_
/* Include Serenity Internal Headers */

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Atom;
/**
 * @class BasisExtension BasisExtension.h
 * @brief A class to systematically extend the basis of subsystems in a sDFT calculation.
 *
 * Currently the total overlap of all active basis functions with the basis function in question
 * is used to determine whether the basis function \f$ \chi_j^B \f$ is added to the basis of the active system:\n\n
 *
 * \f$ \tau^\mathrm{Basis} < \sum_{i\in A^\mathrm{core}}|\langle \chi_i^A|\chi_j^B\rangle|\f$,\n\n
 *
 * where \f$ \tau^\mathrm{Basis}\f$ is the overlap threshold and the sum runs over all all basis functions centered in
 * the active subsystem.
 */
class BasisExtension {
 public:
  /**
   * @brief Default constructor.
   */
  BasisExtension() = default;
  /**
   * @brief Default destructor.
   */
  virtual ~BasisExtension() = default;
  /**
   * @brief Extends the basis of each system via an overlap scheme within the basis of all other systems.
   * @param systems The systems.
   * @param overlapThreshold The overlap threshold.
   */
  void extendAllSystems(std::vector<std::shared_ptr<SystemController>> systems, double overlapThreshold);

 private:
  /// @brief The sets of atoms + ghost atoms for each system.
  std::vector<std::vector<std::shared_ptr<Atom>>> _extendedAtomSets;
  ///@brief The important shells per system and atom.
  std::vector<std::vector<Eigen::VectorXi>> _importantShells;

  /* Helper Functions */
  /**
   * @brief Adds ghost atoms with the additional basis functions to the vector activeAtoms
   *
   * The basis functions centered on environment atoms are checked whether their basis
   * functions may be important for the active system. The basis functions which are
   * important are added to the active system basis by adding the corresponding ghost atoms.
   *
   * @param newImportantShells The important shells, which should be added
   * @param environmentSystem The system controller of the environment system
   * @param activeAtoms The vector of the active system's atoms to which the ghost atoms will be added.
   */
  inline void addGhostAtoms(std::vector<Eigen::VectorXi> newImportantShells, std::shared_ptr<SystemController> environmentSystem,
                            std::vector<Eigen::VectorXi>& activeImportantShells,
                            std::vector<std::shared_ptr<Atom>>& activeAtoms);
  /**
   * @brief Calculates the list of important shells from a given environment system for the active system.
   * @param environmentSystem
   * @param activeSystem
   * @param overlapThreshold The overlap threshold which is used to determine the "importance".
   * @return A vector containing 0 or 1 for each shell of the environment system.
   */
  inline std::vector<Eigen::VectorXi> getListOfImportantShells(std::shared_ptr<SystemController> environmentSystem,
                                                               std::shared_ptr<SystemController> activeSystem,
                                                               double overlapThreshold);

  /**
   * @brief Builds a vector which contains the atoms of the active system.
   * @param activeSystem
   * @return The vector.
   */
  inline std::vector<std::shared_ptr<Atom>> buildActiveSystemAtoms(std::shared_ptr<SystemController> activeSystem);

  /**
   * @brief Builds the basis from the retained ghost atoms+active atoms and their shells.
   * @param activeSystem
   * @param activeAtoms The active atoms + all ghost atoms which are used to expand the basis of the system.
   */
  inline void buildBasisFromGeometry(std::shared_ptr<SystemController> activeSystem,
                                     std::vector<std::shared_ptr<Atom>> activeAtoms,
                                     std::vector<Eigen::VectorXi> importantActiveShells);
};

} /* namespace Serenity */

#endif /* BASIS_BASISEXPANSION_H_ */
