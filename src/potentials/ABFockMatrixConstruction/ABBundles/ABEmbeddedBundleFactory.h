/**
 * @file ABEmbeddedBundleFactory.h
 *
 * @date 18 Aug 2019
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

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLEFACTORY_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLEFACTORY_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"                          //Definition of a DensityMatrixController.
#include "misc/RememberingFactory.h"                                        //Definition of a RememberingFactory.
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundle.h" //The product.
#include "settings/EmbeddingSettings.h"                                     //Definition of embedding settings.
#include "settings/Options.h"                                               //Spin polarization.
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector
namespace Serenity {

/* Forward Declarations */
class SystemController;
class BasisController;
class Geometry;

/**
 * @class ABEmbeddedBundleFactory ABEmbeddedBundleFactory.h
 * @brief A factory that constructs embedded AB potential bundles.\n
 *        This class is a singleton that is constructed on the first call of the
 *        underlying produce function.
 */
template<Options::SCF_MODES SCFMode>
class ABEmbeddedBundleFactory
  : public RememberingFactory<ABEmbeddedBundle<SCFMode>, std::shared_ptr<SystemController>, std::shared_ptr<BasisController>,
                              std::shared_ptr<Geometry>, std::vector<std::shared_ptr<SystemController>>,
                              const std::shared_ptr<EmbeddingSettings>, bool> {
 private:
  /**
   * @brief Construction is handeled indirectly by the produce function.
   */
  ABEmbeddedBundleFactory() = default;

 public:
  /**
   * @brief Gets or produces an ABEmbeddedBundle.
   * @param activeSystem The active system.
   * @param basisControllerB The basis controller defining basis set B.
   * @param geometryB The associated geometry to basis set B.
   * @param environmentSystems The environment systems.
   * @param embeddingSettings The embedding settings.
   * @param topDown Flag for top-down type calculations.
   * @return The AB bundle.
   */
  static std::shared_ptr<ABEmbeddedBundle<SCFMode>>
  produce(std::shared_ptr<SystemController> activeSystem, std::shared_ptr<BasisController> basisControllerB,
          std::shared_ptr<Geometry> geometryB, std::vector<std::shared_ptr<SystemController>> environmentSystems,
          const std::shared_ptr<EmbeddingSettings> embeddingSettings, bool topDown);
  /**
   * @brief Default destructor.
   */
  ~ABEmbeddedBundleFactory() = default;

 private:
  /**
   * @brief Produce a new ABEmbeddedBundle
   * @param activeSystem The active system.
   * @param basisControllerB The basis controller defining basis set B.
   * @param geometryB The associated geometry to basis set B.
   * @param environmentSystems The environment systems.
   * @param embeddingSettings The embedding settings.
   * @param topDown Flag for top-down type calculations.
   * @return The new ABEmbeddedBundle.
   */
  std::unique_ptr<ABEmbeddedBundle<SCFMode>>
  produceNew(std::shared_ptr<SystemController> activeSystem, std::shared_ptr<BasisController> basisControllerB,
             std::shared_ptr<Geometry> geometryB, std::vector<std::shared_ptr<SystemController>> environmentSystems,
             const std::shared_ptr<EmbeddingSettings> embeddingSettings, bool topDown);
  /**
   * @brief Split the environment densities for hybrid approaches.
   * @param activeSystem The active system.
   * @param environmentSystems The environment systems.
   * @param basisFunctionRatio The retained basis function ratio (@see EmbeddingSettings.h).
   * @param borderAtomThreshold The retained border atom threshold (@see EmbeddingSettings.h).
   * @return
   */
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>
  getNotProjectedEnvironmentDensityMatrixControllers(std::shared_ptr<SystemController> activeSystem,
                                                     std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                     double basisFunctionRatio, double borderAtomThreshold);
  /**
   * @brief Construct an AB aux basis controller.
   * @param activeSystem The active system
   * @param geometryB The geometry associated to B.
   * @return The AB aux basis controller.
   */
  std::shared_ptr<BasisController> getABAuxBasisController(std::shared_ptr<SystemController> activeSystem,
                                                           std::shared_ptr<Geometry> geometryB);
};

/**
 * @brief The factories.
 */
static std::unique_ptr<ABEmbeddedBundleFactory<RESTRICTED>> _restrictedFactory;
static std::unique_ptr<ABEmbeddedBundleFactory<UNRESTRICTED>> _unrestrictedFactory;

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABBUNDLES_ABEMBEDDEDBUNDLEFACTORY_H_ */
