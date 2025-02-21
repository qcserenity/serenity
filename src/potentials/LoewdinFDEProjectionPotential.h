/**
 * @file LoewdinFDEProjectionPotential.h
 *
 * @date April 22, 2024
 * @author Denis G. Artiukhin
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

#ifndef POTENTIALS_LOEWDINFDEPROJECTIONPOTENTIAL_H_
#define POTENTIALS_LOEWDINFDEPROJECTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "potentials/Potential.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/Options.h"

namespace Serenity {

/* Forward declarations */
class SystemController;
struct EmbeddingSettings;

/**
 * @class LoewdinFDEProjectionPotential LoewdinFDEProjectionPotential.h
 * @brief A class implementing orbital-dependent Loewdin-type corrections for the non-additive kinetic energy.\n
 *
 *        The idea is based on the use of Slater determinants composed of non-orthogonal molecular orbitals,
 *        the evaluation of kinetic energy expectation values, and the Neumann expansion of the inverse
 * molecular-orbital overlap matrix.\n\n
 *
 * Implementation details and notations follow:\n
 *  J. Chem. Phys. 162, 054117 (2025)\n
 */
template<Options::SCF_MODES SCFMode>
class LoewdinFDEProjectionPotential : public Potential<SCFMode>,
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
   */
  LoewdinFDEProjectionPotential(std::shared_ptr<SystemController> activeSystem,
                                std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                const EmbeddingSettings& settings);
  /**
   * @brief Default destructor.
   */
  virtual ~LoewdinFDEProjectionPotential() = default;

  /**
   * @return The fock matrix for the embedded/active system.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @return Dummy(0) matrix.
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Calculates the associated energy. As opposed to projection-based embedding,
   * this contribution is not equal to zero.
   * @param P The density matrix.
   * @return  The energy value.
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
  /// @brief The potential in matrix representation.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  /// @brief The active system.
  std::weak_ptr<SystemController> _activeSystem;
  /// @brief The environment systems.
  std::vector<std::weak_ptr<SystemController>> _environmentSystems;
  /// @brief The density matrix controllers of the environment systems.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDensityCont;
  // @brief Truncation order for the Loewdin series
  unsigned int _loewdinOrder;
  // @brief Weights (scaling factors) for terms of the Loewdin series
  std::vector<double> _loewdinWeights;

  /// @brief The overlap matrix of the active system basis set with all environment basis sets.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _overlapABs;

  /* Storing intermediate matrices for further re-use */
  // For notations see the ESI of the corresponding article
  SPMatrix<SCFMode> _kAB;
  SPMatrix<SCFMode> _lAA;
  SPMatrix<SCFMode> _mAB;
  SPMatrix<SCFMode> _nAA;

  /* Fock-matrix corrections for energy calculations */
  SPMatrix<SCFMode> _t2Corr;
  SPMatrix<SCFMode> _t3Corr;

  /* Functions performing actual computations of expansion terms*/
  void computeFirstOrderTruncSeries();
  void computeSecondOrderTruncSeries();
  void computeThirdOrderTruncSeries();
};

} /* namespace Serenity */

#endif /* POTENTIALS_LOEWDINFDEPROJECTIONPOTENTIAL_H_ */
