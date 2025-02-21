/**
 * @file HuzinagaProjectionPotential.h
 *
 * @date Nov 23, 2017
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

#ifndef POTENTIALS_HUZINAGAPROJECTIONPOTENTIAL_H_
#define POTENTIALS_HUZINAGAPROJECTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "grid/GridController.h"
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundle.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/Potential.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/SparseCore>

namespace Serenity {

/* Forward declarations */
class SystemController;
struct EmbeddingSettings;

/**
 * @class HuzinagaProjectionPotential HuzinagaProjectionPotential.h
 * @brief A class for handling the Huzinaga operator and Hoffmann's external orthogonality
 *        approach in any FDE like calculation.\n
 *
 *        The operator effectively mirrors the energy levels of the occupied environment orbitals
 *        at 0. The level at which the levels are mirrored can be adjusted by the fermi-shifted
 *        operator.\n\n
 *
 * The Huzinaga operator according to:\n
 *  J. Chem. Theory Comput. 13, 1503-1508 (2017)\n
 *  J. Chem. Phys. 55, 5543 (1971)\n
 *  J. Chem. Phys. 145, 064107 (2016)\n\n
 *
 * Fermi-shifted Huzinaga operator:\n
 *   J. Chem. Theory Comput. 14, 1928-1942 (2018)\n
 */
template<Options::SCF_MODES SCFMode>
class HuzinagaProjectionPotential : public Potential<SCFMode>,
                                    public ObjectSensitiveClass<Basis>,
                                    public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   *
   * This constructor prepares the calculations by building pairs of environment systems with the active system.
   * These pairs are needed in order to calculate the Fock matrix and overlap matrix of the joint system.
   * The constructor expects a guess electronic structure of the active system with SPIN=SCFMode!
   *
   * @param activeSystem The active system.
   * @param environmentSystems The environment systems.
   * @param settings The embedding settings.
   * @param activeFockMatrix The AA potential bundle of the active system.
   * @param topDown Flag for top-down calculations.
   * @param supersystemgrid The supersystem grid controller.
   * @param gridCutOff Optional grid cut off for hybrid approaches.
   * @param allEConts Optional energy component controllers for hybrid approaches.
   * @param fermiShift Optional fermi shift for the operator..
   */
  HuzinagaProjectionPotential(std::shared_ptr<SystemController> activeSystem,
                              std::vector<std::shared_ptr<SystemController>> environmentSystems,
                              const EmbeddingSettings& settings,
                              std::shared_ptr<PotentialBundle<SCFMode>> activeFockMatrix = nullptr, bool topDown = false,
                              std::shared_ptr<GridController> supersystemgrid = nullptr, double gridCutOff = 0.0,
                              std::vector<std::shared_ptr<EnergyComponentController>> allEConts = {nullptr},
                              double fermiShift = 0.0);
  /**
   * @brief Default destructor.
   */
  virtual ~HuzinagaProjectionPotential() = default;

  /**
   * @return The fock matrix for the embedded/active system.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @return Dummy(0) matrix.
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Calculates the associated energy. Only different from 0 if a hybrid approach is used.
   * @param P The density matrix.
   * @return 0. There is no energy associated with this potential.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Potential is linked to the grid and density.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  }

 private:
  ///@brief Sorting/Projection matrices from basis A to basis set B.
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> _AtoBProjections;
  ///@brief Sorting/Projection matrices from the differential basis between B and A (diff  = {B}/{A}) to basis set B.
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> _difftoBProjections;
  /// @brief The potential in matrix representation.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  /// @brief The active system.
  std::weak_ptr<SystemController> _activeSystem;
  /// @brief The environment systems.
  std::vector<std::weak_ptr<SystemController>> _environmentSystems;
  /// @brief The block of the fock matrix from which missing fock matrix elements are extracted if possible.
  std::shared_ptr<PotentialBundle<SCFMode>> _activeFockMatrix;
  /// @brief The density matrix controllers of the environment systems.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDensityCont;
  /// @brief AB-Potentials
  std::vector<std::shared_ptr<ABEmbeddedBundle<SCFMode>>> _abEmbeddedBundles;
  double _fermiShift;
  /*
   * Projection truncation
   */
  /// @brief The overlap matrix of the active system basis set with all environment basis sets.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _s_ABs;
  /// @brief The non additive kinetic energy potential for not projected subsystems (optional).
  std::shared_ptr<NAddFuncPotential<SCFMode>> _naddKinPot;
  /// @brief The density controllers of the not-projected densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _notProjectedEnvDensities;
  /* Helper Functions */
  /**
   * @brief builds the outer diagonal Fock matrix of act. + env[iEnv]
   * @param iEnv The index of the environment system.
   */
  SPMatrix<SCFMode> buildOuterDiagonalFockMatrix(unsigned int iEnv);
  /**
   * @brief Writes the current average overlap of all occupied environmental orbitals with the occ. active system
   * orbitals.
   */
  void writeInterSubsystemOccOverlap();
  /**
   * @brief Adjust the settings for the call to the ABPotentialBundle factory.
   * @param settings The settings to be adjusted.
   * @return A pointer on the manipulated settings.
   */
  std::shared_ptr<EmbeddingSettings> adjustSettings(EmbeddingSettings& settings);
};

} /* namespace Serenity */

#endif /* POTENTIALS_HUZINAGAPROJECTIONPOTENTIAL_H_ */
