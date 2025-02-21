/**
 * @file FDEPotentialBundleFactory.h
 *
 * @author Moritz Bensberg
 * @date Aug 21, 2019
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

#ifndef POTENTIALS_BUNDLES_FDEPOTENTIALBUNDLEFACTORY_H_
#define POTENTIALS_BUNDLES_FDEPOTENTIALBUNDLEFACTORY_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "energies/EnergyComponentController.h"
#include "misc/RememberingFactory.h"
#include "potentials/bundles/PotentialBundle.h"
#include "tasks/FDETask.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class SystemController;
class GridController;
template<Options::SCF_MODES SCFMode>
class PCMPotential;
template<Options::SCF_MODES SCFMode>
class ESIPotentials;
template<Options::SCF_MODES SCFMode>
class ECPInteractionPotential;

template<Options::SCF_MODES SCFMode>
/**
 * @class FDEPotentialBundleFactory FDEPotentialBundleFactory.h
 * @brief Constructs an embedding potential bundle based on a set of EmbeddingSettings and optional
 *        settings.
 */
class FDEPotentialBundleFactory
  : public RememberingFactory<
        PotentialBundle<SCFMode>, std::shared_ptr<SystemController>, std::shared_ptr<DensityMatrixController<SCFMode>>,
        std::vector<std::shared_ptr<SystemController>>, std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>,
        const std::shared_ptr<EmbeddingSettings>, std::shared_ptr<GridController>, std::shared_ptr<SystemController>,
        bool, bool, double, std::vector<std::shared_ptr<EnergyComponentController>>, unsigned int> {
 private:
  /**
   * @brief Default constructor.
   */
  FDEPotentialBundleFactory() = default;

 public:
  /**
   * @brief Default destructor.
   */
  virtual ~FDEPotentialBundleFactory() = default;

  /**
   * @brief Get or produce the potential bundle. Either a FDEPotentials or PBEPotentials object
   *        are generated. The resulting potential bundle will be able to construct the embedded
   *        fock matrix for the active system in the potential of the environment systems.
   * @param activeSystem            The active system controller.
   * @param activeDensMatController The density matrix controller of the active system.
   * @param environmentSystems      The environment systems.
   * @param envDensMatController    The density matrix controller of the environement.
   * @param settings                The EmbeddingSettings.
   * @param grid                    The grid controller.
   * @param supersystem (optional)  The supersystem from a top-down calculation. This is needed for
   *                                top-down potential reconstruction.
   * @param topDown (optional)      Additional flag whether the calculation is top-down or bottom-up.
   * @param noSuperRecon (optional) If true, no potential reconstruction for supersystem is performed.
   *                                Applies to potential reconstruction only .
   * @param gridCutOff (optional)   A flag whether the integration grid is cut off.
   *                                Applies to bottom-up calculations only.
   * @param eConts (optional)       The energy component controller of the environment systems.
   *                                Applies to bottom-up calculations only.
   * @return The potential bundle. If top-down=true -> PBEPotentials, else->FDEPotentials.
   */
  static std::shared_ptr<PotentialBundle<SCFMode>>
  produce(std::shared_ptr<SystemController> activeSystem,
          std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
          std::vector<std::shared_ptr<SystemController>> environmentSystems,
          std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
          const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
          std::shared_ptr<SystemController> supersystem = nullptr, bool topDown = false, bool noSuperRecon = true,
          double gridCutOff = -1.0, std::vector<std::shared_ptr<EnergyComponentController>> eConts = {},
          unsigned int firstPassiveSystemIndex = 99999);

 private:
  // The actual worker function. For documentaion see produce(..)
  std::unique_ptr<PotentialBundle<SCFMode>>
  produceNew(std::shared_ptr<SystemController> activeSystem,
             std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
             std::vector<std::shared_ptr<SystemController>> environmentSystems,
             std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
             const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
             std::shared_ptr<SystemController> supersystem, bool topDown, bool noSuperRecon, double gridCutOff,
             std::vector<std::shared_ptr<EnergyComponentController>> eConts, unsigned int firstPassiveSystemIndex);
  // Building the mixed embedding potential for exact/approx embedding
  static std::unique_ptr<PotentialBundle<SCFMode>>
  buildMixedEmbedding(std::shared_ptr<SystemController> activeSystem,
                      std::shared_ptr<DensityMatrixController<SCFMode>> activeDensMatController,
                      std::vector<std::shared_ptr<SystemController>> environmentSystems,
                      std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensMatController,
                      const std::shared_ptr<EmbeddingSettings> settings, std::shared_ptr<GridController> grid,
                      bool topDown, double gridCutOff, std::vector<std::shared_ptr<EnergyComponentController>> eConts,
                      std::shared_ptr<PotentialBundle<SCFMode>> potBundle,
                      std::shared_ptr<PotentialBundle<SCFMode>> esiPot, std::shared_ptr<PCMPotential<SCFMode>> pcm,
                      std::shared_ptr<ECPInteractionPotential<SCFMode>> ecpInt_total);
};

static std::unique_ptr<FDEPotentialBundleFactory<RESTRICTED>> _restrictedFDEFactory;
static std::unique_ptr<FDEPotentialBundleFactory<UNRESTRICTED>> _unrestrictedFDEFactory;

} /* namespace Serenity */

#endif /* POTENTIALS_BUNDLES_FDEPOTENTIALBUNDLEFACTORY_H_ */
