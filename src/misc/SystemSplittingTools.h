/**
 * @file SystemSplittingTools.h
 *
 * @date Sep 25, 2018
 * @author Moritz Bensberg
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

#ifndef MISC_SYSTEMSPLITTINGTOOLS_H_
#define MISC_SYSTEMSPLITTINGTOOLS_H_
/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "data/matrices/SPMatrix.h"
#include "math/Matrix.h"
#include "data/matrices/DensityMatrixController.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>
#include <Eigen/Dense>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Atom;
class Geometry;
template<Options::SCF_MODES T> class ElectronicStructure;
template<Options::SCF_MODES T> class MatrixInBasis;
class BasisController;


/**
 * @class SystemSplittingTools SystemSplittingTools.h
 * @brief A collection of tools for splitting and adding systems or associated properties.
 */
template<Options::SCF_MODES SCFMode>
class SystemSplittingTools {
private:
  /* Purely static. Never instantiated. */
  SystemSplittingTools()=default;
  virtual ~SystemSplittingTools()=default;
public:
  /**
   * @brief Splits the electronic structure of the system based on the selected orbitals for the active system.
   * @param system The system which will be split.
   * @param activeOrbtials The selection of the active orbitals.
   * @return Electronic structures for the subsystems (act : env)
   */
  static std::pair<std::shared_ptr<ElectronicStructure<SCFMode> >,std::shared_ptr<ElectronicStructure<SCFMode> > >
  splitElectronicStructure(
      std::shared_ptr<SystemController> system,
      SpinPolarizedData<SCFMode,std::vector<bool> >& activeOrbtials);

  /**
   * @brief Construct an atom from a dummy.
   * @param atom The dummy-atom
   * @return The non-dummy atom.
   */
  static std::shared_ptr<Atom> makeAtomFromDummy(std::shared_ptr<Atom> atom);

  /**
   * @brief Partitions the occupied orbitals based on their <population> on the <activeAtoms>.
   * @param orbitalPopulations The populations of the orbitals on all atoms.
   * @param supersystem The systems which orbitals will be partitioned.
   * @param activeAtoms The atoms of the system which are considered to be active.
   * @param orbitalThreshold The localization threshold for the orbital picking.
   * @param nOccOrbsActive (optional) The maximum number of orbitals which should be picked.
   *                       (needs enforceCharges=true)
   * @param enforceCharges (optional) Enforce a specific number of orbitals which are picked.
   * @return A list of booleans which are true if the orbital was picked.
   */
  static SpinPolarizedData<SCFMode,std::vector<bool> > partitionOrbitals(
      SPMatrix<SCFMode>& orbitalPopulations,
      std::shared_ptr<SystemController> supersystem,
      std::vector<bool> activeAtoms,
      double orbitalThreshold,
      SpinPolarizedData<SCFMode,unsigned int> nOccOrbsActive = SpinPolarizedData<SCFMode,unsigned int> (0),
      bool enforceCharges =false);

  /**
   * @brief Projects a matrix into a new basis.
   * @param oldMatrix The matrix which should be projected.
   * @param newBasis The new basis.
   * @param newOverlap (optional) The overlap matrix of the new basis.
   * @return The projected matrix.
   */
  static MatrixInBasis<SCFMode> projectMatrixIntoNewBasis(
      const MatrixInBasis<SCFMode>& oldMatrix,
      std::shared_ptr<BasisController> newBasis,
      std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > newOverlap = nullptr);

  /**
   * @brief Calculates the occupations of the subsystems from the orbital selection.
   * @param activeOrbtials The selection of active orbitals.
   * @return The occupation numbers for environment and active system (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> >
  getNOccupiedOrbitals(SpinPolarizedData<SCFMode,std::vector<bool> >& activeOrbtials);

  /**
   * @brief Calculates the number of electrons in active and environment system.
   * @param activeOrbtials List of occupied orbitals of the active system. Length of the vector
   *                       denotes the number of active orbitals and True/False denotes which
   *                       orbital belongs to the active system.
   * @return The number of electrons for the active and environment system. (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> >
  getNElectrons(SpinPolarizedData<SCFMode,std::vector<bool> >& activeOrbtials) {
    return getNElectrons(getNOccupiedOrbitals(activeOrbtials));
  }
  /**
   * @brief Calculates the number of electrons in active and environment system.
   * @param occ Number of occupied orbitals in the active and environment system.
   * @return The number of electrons for the active and environment system. (act : env)
   */
  static std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> >
  getNElectrons(std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> > occ);

  /**
   * @brief Returns the access of alpha orbitals.
   * @param nOcc The occupations.
   * @return The spin corresponding to the occupations.
   */
  static int getSpin(SpinPolarizedData<SCFMode,unsigned int> nOcc);

  /**
   * @brief Get the index of the atom in the given geometry. The dummy attribute is ignored for this
   *         matching.
   * @param geometry The geometry.
   * @param atom The environment atom.
   * @return The index of the atom in the given geometry. If the atom can not be matched, the returned index
   *         will be larger than the number of atoms in the geometry.
   */
  static unsigned int matchAtom(
      std::shared_ptr<Geometry> geometry,
      std::shared_ptr<Atom> atom);

  /**
   * @brief Selects distant orbitals based on the active system basis set.
   * @param orbitalPopulations The atom-wise orbital populations.
   * @param activeSystem The active system.
   * @param environmentSystem The environment system in question.
   * @param basisFunctionRatio The basis function shell ratio for the selection of distant atoms.
   * @param borderAtomThreshold The threshold for the selection of the distant orbitals.
   * @return Individual flags for the distant orbitals.
   */
  static SpinPolarizedData<SCFMode,std::vector<bool> >selectDistantOrbitals(
      SPMatrix<SCFMode>& orbitalPopulations,
      std::shared_ptr<SystemController> activeSystem,
      std::shared_ptr<SystemController> environmentSystem,
      double basisFunctionRatio,
      double borderAtomThreshold);

  /**
   * @brief Construct the density matrix of the orbitals not included in the projection operator for
   *        a specific system.
   * @param environmentSystem The environment system.
   * @param distantOrbitals The selection of distant orbitals.
   * @return The corresponding density matrix.
   */
  static std::shared_ptr<DensityMatrix<SCFMode> > buildNonOrthogonalDensityMatrix(
      std::shared_ptr<SystemController> environmentSystem,
      SpinPolarizedData<SCFMode,std::vector<bool> > distantOrbitals);

  /**
   * @brief Calculates the overlap between the basis sets of active and environment systems.
   *        A nullptr is returned if the total overlap is below the given threshold.
   * @param activeSystem The active systems.
   * @param environmentSystem The environment systems.
   * @param truncationThreshold Total overlap threshold.
   * @return The overlap matrices or nullptr in case of a not projected system.
   */
  static std::vector<std::shared_ptr<Eigen::MatrixXd> > getProjectedSubsystems(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems,
      double truncationThreshold = -10);

  /**
   * @brief Constructs the density matrix controllers with the given SCFMode.
   * @param environmentSystems The environment system controllers.
   * @param topDown A flag for top-down calculations.
   * @return The associated density matrix controllers with the correct SCFMode.
   */
  static std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > getEnvironmentDensityControllers(
      std::vector<std::shared_ptr<SystemController> > environmentSystems,
      bool topDown = false);
};

} /* namespace Serenity */

#endif /* MISC_SYSTEMSPLITTINGTOOLS_H_ */
