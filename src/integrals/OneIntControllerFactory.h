/**
 * @file   OneIntControllerFactory.h
 *
 * @date   last rework Nov 29. 2016
 * @author Thomas Dresselhaus, Jan Unsleber
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
#ifndef ONEINTCONTROLLERFACTORY_H
#define ONEINTCONTROLLERFACTORY_H
/* Include Serenity Internal Headers */
#include "integrals/OneElectronIntegralController.h"
#include "misc/RememberingFactory.h"

namespace Serenity {
/* Forward declarations */
/**
 * @class OneIntControllerFactory OneIntControllerFactory.h
 *
 * @brief This factory produces OneElectronIntegralControllers with a given basis
 *        and geometry.
 *
 * To not calculate redundant data pointers to all the created instances are stored and thus no
 * second controller will be produced for the same geometry and basis.
 * This class is a singleton
 */
class OneIntControllerFactory
  : public RememberingFactory<OneElectronIntegralController, const std::shared_ptr<BasisController>,
                              const std::shared_ptr<const Geometry>, const std::shared_ptr<ExternalChargeController>> {
 public:
  /**
   * @brief One of two singleton 'Constructors'
   *
   * @return A reference to the only existing version of the OneIntControllerFactory object.
   */
  static OneIntControllerFactory& getInstance() {
    return *OneIntControllerFactory::getSharedPtr();
  }
  /**
   * @brief One of two Singleton 'Constructors'
   * @return A shared pointer to the only existing version of the OneIntControllerFactory object.
   */
  static std::shared_ptr<OneIntControllerFactory> getSharedPtr() {
    static std::shared_ptr<OneIntControllerFactory> instance(new OneIntControllerFactory);
    return instance;
  }
  /**
   * @brief Deleted -> singleton.
   */
  OneIntControllerFactory(OneIntControllerFactory const&) = delete;
  /**
   * @brief Deleted -> singleton.
   */
  void operator=(OneIntControllerFactory const&) = delete;
  /**
   * @brief Default destructor.
   */
  virtual ~OneIntControllerFactory() = default;
  /**
   * @brief Controllers are precached
   *
   * , i.e. only with a new basis an actually new controller is produced. In any case the garbage
   * collector is already notified that the calling class uses the returned instance.
   *
   * @param basisController    The basis controller defining the AO basis.
   * @param geometry           The molecular geometry, necessary for nuclei-electron attraction integrals.
   * @param externalCharges    The controller holding the external charges (if any).
   * @returns a new instance or the pointer to an old instance if it already exists.
   */
  std::shared_ptr<OneElectronIntegralController> produce(std::shared_ptr<BasisController> basisController,
                                                         std::shared_ptr<const Geometry> geometry,
                                                         std::shared_ptr<ExternalChargeController> externalCharges) {
    return getOrProduce(basisController, geometry, externalCharges);
  }

 private:
  /**
   * @brief Private destructor -> singleton.
   */
  OneIntControllerFactory() = default;

  std::unique_ptr<OneElectronIntegralController>
  produceNew(const std::shared_ptr<BasisController> basisController, const std::shared_ptr<const Geometry> geometry,
             const std::shared_ptr<ExternalChargeController> externalCharges) override final {
    return std::unique_ptr<OneElectronIntegralController>(
        new OneElectronIntegralController(basisController, geometry, externalCharges));
  }
};

} /* namespace Serenity */
#endif /* ONEINTCONTROLLERFACTORY_H */
