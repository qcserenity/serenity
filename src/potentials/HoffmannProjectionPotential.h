/**
 * @file HoffmannProjectionPotential.h
 *
 * @author Moritz Bensberg
 * @date Aug 26, 2019
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

#ifndef POTENTIALS_HOFFMANNPROJECTIONPOTENTIAL_H_
#define POTENTIALS_HOFFMANNPROJECTIONPOTENTIAL_H_
/* Include Serenity Internal Headers */
#include "basis/Basis.h"                           //ObjectSensitiveClass
#include "data/matrices/DensityMatrixController.h" //DensityMatrixController definition
#include "potentials/Potential.h"                  //Potential definition.
#include "potentials/bundles/PotentialBundle.h"    //Diagonal block of the fock matrix.
#include "settings/EmbeddingSettings.h"            //EmbeddingSettings definition.
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {
/* Forward declarations */
class SystemController;
class EnergyComponentController;
class GridController;

/**
 * @class HoffmannProjectionPotential HoffmannProjectionPotential.h
 * @brief A class for the calculation of a Hoffmann-type projection operator.\n\n
 *
 * Hoffmann's external orthogonality approach:\n
 * According to:\n
 *  Theory. Annu. Rep. Comput. Chem. 8,53-70 (2012) and \n
 *  J. Phys. Chem. A 118, 9182–9200 (2014).\n
 */
template<Options::SCF_MODES SCFMode>
class HoffmannProjectionPotential : public Potential<SCFMode>,
                                    public ObjectSensitiveClass<Basis>,
                                    public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem         The active system.
   * @param environmentSystems   The environment systems.
   * @param settings             The embedding settings.
   * @param activeFockMatrix     The active fock matrix block without projection as a bundle.
   * @param topDown              Flag for top-down calculations.
   * @param supersystemgrid      The supersystem grid controller.
   * @param gridCutOff           Optional grid cut off for hybrid approaches.
   * @param allEConts            Optional energy component controller for hybrid approaches.
   */
  HoffmannProjectionPotential(std::shared_ptr<SystemController> activeSystem,
                              std::vector<std::shared_ptr<SystemController>> environmentSystems,
                              const EmbeddingSettings& settings,
                              std::shared_ptr<PotentialBundle<SCFMode>> activeFockMatrix = nullptr, bool topDown = false,
                              std::shared_ptr<GridController> supersystemgrid = nullptr, double gridCutOff = 0.0,
                              std::vector<std::shared_ptr<EnergyComponentController>> allEConts = {});
  /**
   * @brief Defualt destrutor.
   */
  virtual ~HoffmannProjectionPotential() = default;

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
  /**
   * @brief Adjust the settings for the call to the ABPotentialBundle factory.
   * @param settings The settings to be adjusted.
   * @return A pointer on the manipulated settings.
   */
  std::shared_ptr<EmbeddingSettings> adjustSettings(EmbeddingSettings& settings);
  /**
   * @brief Get the non picked system controllers.
   * @param pickedSystem The picked system controller.
   * @return The non picked system controllers.
   */
  std::vector<std::shared_ptr<SystemController>> getOtherSystemController(std::shared_ptr<SystemController> pickedSystem);
  /**
   *
   * @brief Get the non picked density matrix controllers.
   * @param pickedDensityMatrixController
   * @param activeDensityMatrixController
   * @param environmentDensityMatrixController
   * @return The non picked density matrix controllers.
   */
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>
  getOtherDensityMatrixController(std::shared_ptr<DensityMatrixController<SCFMode>> pickedDensityMatrixController,
                                  std::shared_ptr<DensityMatrixController<SCFMode>> activeDensityMatrixController,
                                  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> environmentDensityMatrixController);

  /*
   * @brief Split the construction of this potential so that the FDEPotentialBundleFactory can be a RememberingFactory.
   * Since getOrProduce() is locked when this Potential is built, its own calls to the Factory are now done in this
   * function after construction. It is invoked the first time getMatrix() is called.
   */
  void finishSetup();

  /**
   * @brief The Huzinaga projection part.
   */
  std::shared_ptr<Potential<SCFMode>> _huzinagaProjection;
  /**
   * @brief The diagonal fock matrix blocks as bundles.
   */
  std::vector<std::shared_ptr<PotentialBundle<SCFMode>>> _embeddingBundle;
  /**
   * @brief The active system controller.
   */
  std::weak_ptr<SystemController> _activeSystem;
  /**
   * @brief The environment system controllers.
   */
  std::vector<std::weak_ptr<SystemController>> _environmentSystems;
  /// @brief The overlap matrix of the active system basis set with all environment basis sets.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _s_ABs;
  /// @brief The density controllers of the not-projected densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _notProjectedEnvDensities;
  /// @brief The density matrix controllers of the environment systems.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _envDensityCont;
  /// @brief The potential in matrix representation.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  /// @brief The EmbeddingSettings after the adjustSettings() function has been applied.
  std::shared_ptr<EmbeddingSettings> _adjustedSettings_ptr;
  /// @brief Bool whether the calculation is top-down or bottom-up.
  bool _topDown;
  /// @brief A GridController for the Supersystem integration grid, if it is employed. Otherwise this is a nullptr.
  std::shared_ptr<GridController> _supersystemGrid;
  /// @brief Bool to make sure the post-construction function finishSetup() is only called once.
  bool _finishedSetup = false;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HOFFMANNPROJECTIONPOTENTIAL_H_ */
