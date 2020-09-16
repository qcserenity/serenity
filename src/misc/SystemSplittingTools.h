/**
 * @file SystemSplittingTools.h
 *
 * @date Sep 25, 2018
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

#ifndef MISC_SYSTEMSPLITTINGTOOLS_H_
#define MISC_SYSTEMSPLITTINGTOOLS_H_
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "data/matrices/SPMatrix.h"
#include "math/Matrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Atom;
class Geometry;
template<Options::SCF_MODES T>
class ElectronicStructure;
template<Options::SCF_MODES T>
class MatrixInBasis;
class BasisController;
class OneElectronIntegralController;

/**
 * @class SystemSplittingTools SystemSplittingTools.h
 * @brief A collection of tools for splitting and adding systems or associated properties.
 */
template<Options::SCF_MODES SCFMode>
class SystemSplittingTools {
 private:
  /* Purely static. Never instantiated. */
  SystemSplittingTools() = default;
  virtual ~SystemSplittingTools() = default;

 public:
  /**
   * @brief Splits the electronic structure of the system based on the selected orbitals for the active system.
   * @param system The system which will be split.
   * @param activeOrbtials The selection of the active orbitals.
   * @return Electronic structures for the subsystems (act : env)
   */
  static std::pair<std::shared_ptr<ElectronicStructure<SCFMode>>, std::shared_ptr<ElectronicStructure<SCFMode>>>
  splitElectronicStructure(std::shared_ptr<SystemController> system,
                           SpinPolarizedData<SCFMode, std::vector<bool>>& activeOrbtials);
  /**
   * @brief Resorts the orbital coefficients from the old basis set to the new basis set. Note that basis functions
   *        which are not present in both basis sets will be omitted/lead to zero entries.
   * @param electronicStructure The "old" electronic structure.
   * @param newBasisController The new basis set to express the orbitals in.
   * @return The "new" electronic structure in which the orbital coefficients have been resorted to be expressed in the
   *         new basis set.
   */
  static std::shared_ptr<ElectronicStructure<SCFMode>>
  resortBasisSetOfElectronicStructure(std::shared_ptr<ElectronicStructure<SCFMode>> electronicStructure,
                                      std::shared_ptr<BasisController> newBasisController,
                                      std::shared_ptr<OneElectronIntegralController> oneElectronIntegralController);
  /**
   * @brief Splits the geometry of the system up based on the localization of the (split) orbitals.
   *        Atoms with a total population of active orbitals larger than the given threshold are
   *        assigned to the active system.
   * @param system The system which will be split.
   * @param activeES The electronic structure of the selected environment system.
   * @param localizationThreshold The threshold for the atom assignment.
   * @return Geometries of the active and environment system (act : env)
   */
  static std::pair<std::shared_ptr<Geometry>, std::shared_ptr<Geometry>>
  splitGeometry(std::shared_ptr<SystemController> system, std::shared_ptr<ElectronicStructure<SCFMode>> activeES,
                double localizationThreshold);

  /// @brief Produce a new dummy atom from a non-ghost atom.
  /// @param atom The atom.
  /// @return The new dummy atom.
  static std::shared_ptr<Atom> makeDummyAtom(std::shared_ptr<Atom> atom);
  /**
   * @brief Produce a new atom from a dummy atom.
   * @param atom The atom.
   * @return The non-dummy atom.
   */
  static std::shared_ptr<Atom> makeAtomFromDummy(std::shared_ptr<Atom> atom);

  /**
   * @brief Projects a matrix into a new basis.
   * @param oldMatrix The matrix which should be projected.
   * @param newBasis The new basis.
   * @param newOverlap (optional) The overlap matrix of the new basis.
   * @return The projected matrix.
   */
  static MatrixInBasis<SCFMode>
  projectMatrixIntoNewBasis(const MatrixInBasis<SCFMode>& oldMatrix, std::shared_ptr<BasisController> newBasis,
                            std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED>> newOverlap = nullptr);

  /**
   * @brief Calculates the occupations of the subsystems from the orbital selection.
   * @param activeOrbtials The selection of active orbitals.
   * @return The occupation numbers for environment and active system (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
  getNOccupiedOrbitals(SpinPolarizedData<SCFMode, std::vector<bool>>& activeOrbitals);

  /**
   * @brief Calculates the number of electrons in active and environment system.
   * @param activeOrbtials List of occupied orbitals of the active system. Length of the vector
   *                       denotes the number of active orbitals and True/False denotes which
   *                       orbital belongs to the active system.
   * @return The number of electrons for the active and environment system. (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
  getNElectrons(SpinPolarizedData<SCFMode, std::vector<bool>>& activeOrbtials) {
    return getNElectrons(getNOccupiedOrbitals(activeOrbtials));
  }
  /**
   * @brief Calculates the number of electrons in active and environment system.
   * @param occ Number of occupied orbitals in the active and environment system.
   * @return The number of electrons for the active and environment system. (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
  getNElectrons(std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>> occ);

  /**
   * @brief Returns the access of alpha orbitals.
   * @param nOcc The occupations.
   * @return The spin corresponding to the occupations.
   */
  static int getSpin(SpinPolarizedData<SCFMode, unsigned int> nOcc);

  /**
   * @brief Get the index of the atom in the given geometry. The dummy attribute is ignored for this
   *         matching.
   * @param geometry The geometry.
   * @param atom The environment atom.
   * @return The index of the atom in the given geometry. If the atom can not be matched, the returned index
   *         will be larger than the number of atoms in the geometry.
   */
  static unsigned int matchAtom(std::shared_ptr<Geometry> geometry, std::shared_ptr<Atom> atom);

  /**
   * @brief Selects distant orbitals based on the active system basis set.
   * @param orbitalPopulations The atom-wise orbital populations.
   * @param activeSystem The active system.
   * @param environmentSystem The environment system in question.
   * @param basisFunctionRatio The basis function shell ratio for the selection of distant atoms.
   * @param borderAtomThreshold The threshold for the selection of the distant orbitals.
   * @return Individual flags for the distant orbitals.
   */
  static SpinPolarizedData<SCFMode, std::vector<bool>>
  selectDistantOrbitals(SPMatrix<SCFMode>& orbitalPopulations, std::shared_ptr<SystemController> activeSystem,
                        std::shared_ptr<SystemController> environmentSystem, double basisFunctionRatio,
                        double borderAtomThreshold);

  /**
   * @brief Construct the density matrix of the orbitals not included in the projection operator for
   *        a specific system.
   * @param environmentSystem The environment system.
   * @param distantOrbitals The selection of distant orbitals.
   * @return The corresponding density matrix.
   */
  static std::shared_ptr<DensityMatrix<SCFMode>>
  buildNonOrthogonalDensityMatrix(std::shared_ptr<SystemController> environmentSystem,
                                  SpinPolarizedData<SCFMode, std::vector<bool>> distantOrbitals);
  /**
   * @brief Get the matrix block which is expanded over the basis functions centered on the atoms
   *        with indices atomsA x atomsB.
   * @param matrix The matrix from which the block is extracted.
   * @param atomsA The "row" atom indices.
   * @param atomsB The "column" atom indices.
   * @param atomBasisIndices The map between shells and atoms.
   * @return The matrix block.
   */
  static SPMatrix<SCFMode> getMatrixBlock(const SPMatrix<SCFMode>& matrix, const std::vector<unsigned int> atomsA,
                                          const std::vector<unsigned int> atomsB,
                                          const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices);
  /**
   * @brief Get the matrix block which is expanded over the shells which are
   *        not zero in the given sparse vectors.
   * @param matrix The matrix from which the block is extracted.
   * @param shellsA The "row" shell indices.
   * @param shellsB The "column" shell indices.
   * @return The matrix block.
   */
  static SPMatrix<SCFMode> getMatrixBlockShellWise(const MatrixInBasis<SCFMode>& matrix, const Eigen::SparseVector<int> shellsA,
                                                   const Eigen::SparseVector<int> shellsB);
  /**
   * @brief Diagonalize a matrix in a linear-dependent basis.
   *
   *      The idea to diagonalize a matrix in a basis that is expanded in the AO basis,
   *      but does not span it fully and may be linear dependent.
   *
   * @param R_ij The AO coefficients of the basis.
   * @param s_AO The AO overlap matrix.
   * @param f_AO The matrix in AO basis.
   * @param paoOrthogonalizationThreshold The threshold for the canonical orthogonalization.
   * @param eigenvalues The matrix eigenvalues.
   * @param transformation The transformation to the non-linear dependet basis.
   */
  static void diagonalizationInNonRedundantPAOBasis(const Eigen::MatrixXd& R_ij, const Eigen::MatrixXd& s_AO,
                                                    const Eigen::MatrixXd& f_AO, double paoOrthogonalizationThreshold,
                                                    Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& transformation);
  /**
   * @brief The matrix is squared and checked for values larger than the threshold.
   *        A sparse map is constructed that contains only elements that fullfill this
   *        condition.
   * @param matrix       The matrix.
   * @param mnpThreshold The threshold.
   * @return The sparse map.
   */
  static Eigen::SparseMatrix<int> reduceMatrixToMullikenNetPopulationMap(const Eigen::MatrixXd& matrix, double mnpThreshold);

  /**
   * @brief Calculates the overlap between the basis sets of active and environment systems.
   *        A nullptr is returned if the total overlap is below the given threshold.
   * @param activeSystem The active systems.
   * @param environmentSystem The environment systems.
   * @param truncationThreshold Total overlap threshold.
   * @return The overlap matrices or nullptr in case of a not projected system.
   */
  static std::vector<std::shared_ptr<Eigen::MatrixXd>>
  getProjectedSubsystems(std::shared_ptr<SystemController> activeSystem,
                         std::vector<std::shared_ptr<SystemController>> environmentSystems, double truncationThreshold = -10);

  /**
   * @brief Constructs the density matrix controllers with the given SCFMode.
   * @param environmentSystems The environment system controllers.
   * @param topDown A flag for top-down calculations.
   * @return The associated density matrix controllers with the correct SCFMode.
   */
  static std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>
  getEnvironmentDensityControllers(std::vector<std::shared_ptr<SystemController>> environmentSystems, bool topDown = false);
};

} /* namespace Serenity */

#endif /* MISC_SYSTEMSPLITTINGTOOLS_H_ */
