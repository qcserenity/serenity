/**
 * @file ESIPotentials.h
 *
 * @date Nov 24, 2016
 * @author: Jan Unsleber
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

#ifndef POTENTIALS_BUNDLES_ESIPOTENTIALS_H_
#define POTENTIALS_BUNDLES_ESIPOTENTIALS_H_

/* Include Serenity Internal Headers */

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "geometry/Geometry.h"
#include "potentials/Potential.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/LocalizationOptions.h"
#include "settings/Options.h"

namespace Serenity {
/* Forward declaration */
class SystemController;
/**
 * @class ESIPotentials ESIPotentials.h
 * @brief The electrostatic interaction potentials/energies bundled into one class.
 *
 * This class combines the description of the following interactions between one
 * active system and a multitude of environment systems:
 *  - the nuclear-electron interactions
 *    (both active nucleii with env. electrons and active electrons with env. nucleii)
 *  - the Coulomb repulsion between the active electrons and all environment electrons
 *    (excluding electron-electron repulsion within and between the different
 *     environment systems)
 *  - the Coulomb repulsion between the active system nucleii and the environment nucleii
 *    (excluding nucleus-nucleus repulsion within and between the different
 *     environment systems)
 *  - due to the technical relation of Coulomb and HF-exchange contributions
 *    this class also allows for a consideration of exchange interactions.
 *    (Please see the appropriate constructor for more information.)
 */
template<Options::SCF_MODES SCFMode>
class ESIPotentials : public PotentialBundle<SCFMode>,
                      public ObjectSensitiveClass<Atom>,
                      public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   *
   * This constructor builds a version of the Potential that does not use the RI approximation.
   * The DensityMatrix and Geometry arguments in addition to the systems
   * are required to allow for custom matrices and geometries that are not stored in the
   * systems.
   *
   * This constructor allows for the coupling of the active system to the environment
   * using a HF-like exchange interaction contribution.
   * IMPORTANT: the implementation assumes pairwise orthogonal orbitals between the active
   *            and all environment orbitals. The implementation does NOT solve Loewdin's
   *            general formula for non-orthogonal molecular orbitals.
   *
   * @param actSystem               The single active system all other systems act upon.
   * @param envSystems              The environment systems that are to be coupled to the active system.
   * @param activeDMat              The DensityMatrix of the active system.
   * @param activeGeom              The Geometry of the active system.
   * @param envDMats                The DensityMatrices of the environment systems.
   * @param envGeoms                The geometries of the environment systems.
   * @param firstPassiveSystemIndex First index of a passive environment systems that will not be changed during
   *                                freeze-and-thaw iterations. Used for caching.
   */
  ESIPotentials(std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
                std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat, std::shared_ptr<const Geometry> activeGeom,
                std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                std::vector<std::shared_ptr<const Geometry>> envGeoms, unsigned int firstPassiveSystemIndex = 99999,
                bool useCharges = false,
                Options::POPULATION_ANALYSIS_ALGORITHMS chargeModel = Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN);

  /**
   * @brief Constructor
   *
   * This constructor builds a version of the Potential that uses the RI approximation.
   * The DensityMatrix, auxiliary basis and Geometry arguments in addition to the systems
   * are required to allow for custom matrices etc. that are not stored in the
   * systems.
   *
   * With this constructor it is not possible to add an exact exchange interaction part.
   *
   * @param actSystem               The single active system all other systems act upon.
   * @param envSystems              The environment systems that are to be coupled to the active system.
   * @param activeDMat              The DensityMatrix of the active system.
   * @param activeGeom              The Geometry of the active system.
   * @param envDMats                The DensityMatrices of the environment systems.
   * @param envGeoms                The geometries of the environment systems.
   * @param actAuxBasis             The auxiliary basis set of the active system.
   * @param envAuxBasis             The auxiliary basis sets of the environment systems.
   * @param firstPassiveSystemIndex First index of a passive environment systems that will not be changed during
   *                                freeze-and-thaw iterations. Used for caching.
   */
  ESIPotentials(std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
                std::shared_ptr<DensityMatrixController<SCFMode>> activeDMat, std::shared_ptr<const Geometry> activeGeom,
                std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDMats,
                std::vector<std::shared_ptr<const Geometry>> envGeoms, const std::shared_ptr<BasisController> actAuxBasis,
                std::vector<std::shared_ptr<BasisController>> envAuxBasis, unsigned int firstPassiveSystemIndex = 99999,
                bool useCharges = false,
                Options::POPULATION_ANALYSIS_ALGORITHMS chargeModel = Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN);

  /// @brief Default destructor.
  virtual ~ESIPotentials() = default;

  /**
   * @brief A function to get the entire Fock matrix.
   * @param P The density matrix.
   * @param energies The controller to add all the energy contributions to.
   * @return Returns the current Fock matrix (rebuilt on every call).
   *         (Note that this does not imply that each Potential is rebuilt, they
   *          may very well be cached. But all potentials contributing are added
   *          together in every call)
   */
  FockMatrix<SCFMode> getFockMatrix(const DensityMatrix<SCFMode>& P,
                                    std::shared_ptr<EnergyComponentController> energies) override final;

  /**
   * @brief Getter for the nuclear gradients.
   * @return Returns the combined gradients of all underlying potentials.
   */
  Eigen::MatrixXd getGradients() override final;

  /**
   * @brief The Potential is linked to the environment density and the geometries.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _enAttr = nullptr;
    _nnRep = nullptr;
  };

 private:
  std::shared_ptr<SystemController> _actSystem;
  std::vector<std::shared_ptr<SystemController>> _envSystems;
  std::shared_ptr<DensityMatrixController<SCFMode>> _activeDMat;
  std::shared_ptr<const Geometry> _activeGeom;
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDMats;
  std::vector<std::shared_ptr<const Geometry>> _envGeoms;
  std::shared_ptr<Potential<SCFMode>> _nePot;
  std::shared_ptr<Potential<SCFMode>> _cePot;
  std::shared_ptr<Potential<SCFMode>> _coulPot;
  std::shared_ptr<Potential<SCFMode>> _passiveCoulPot;
  std::unique_ptr<double> _enAttr;
  std::unique_ptr<double> _nnRep;
};

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_ESIPOTENTIALS_H_ */
