/**
 * @file RI_J_IntegralControllerFactory.h
 *
 * @date Mar 9, 2016
 * @author: Kevin Klahr, Jan Unsleber
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

#ifndef INTEGRALS_RI_J_INTEGRALCONTROLLERFACTORY_H_
#define INTEGRALS_RI_J_INTEGRALCONTROLLERFACTORY_H_
/* Include Serenity Internal Headers */
#include "integrals/RI_J_IntegralController.h"
#include "misc/RememberingFactory.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class BasisController;
class ERIControllerFactory;
class Libint;
class MatrixInverter;
class RI_J_IntegralController;

/**
 *  @class RI_J_IntegralControllerFactory RI_J_IntegralControllerFactory.h
 *  @brief This class is a singleton to generate RI_J_IntegralController.
 */
class RI_J_IntegralControllerFactory
  : public RememberingFactory<RI_J_IntegralController, std::shared_ptr<BasisController>, std::shared_ptr<BasisController>,
                              std::shared_ptr<BasisController>, const LIBINT_OPERATOR, const double> {
 public:
  /**
   * @brief One of two singleton 'Constructors'
   *
   * @return A reference to the only existing version of the RI_J_IntegralControllerFactory object.
   */
  static RI_J_IntegralControllerFactory& getInstance() {
    return *RI_J_IntegralControllerFactory::getSharedPtr();
  }
  /**
   * @brief One of two Singleton 'Constructors'
   * @return A shared pointer to the only existing version of the RI_J_IntegralControllerFactory object.
   */
  static std::shared_ptr<RI_J_IntegralControllerFactory> getSharedPtr() {
    static std::shared_ptr<RI_J_IntegralControllerFactory> instance(new RI_J_IntegralControllerFactory);
    return instance;
  }
  /**
   * @brief Deleted -> singleton.
   */
  RI_J_IntegralControllerFactory(RI_J_IntegralControllerFactory const&) = delete;
  /**
   * @brief Deleted -> singleton.
   */
  void operator=(RI_J_IntegralControllerFactory const&) = delete;
  /**
   * @brief Default destructor.
   */
  virtual ~RI_J_IntegralControllerFactory() = default;

  std::shared_ptr<RI_J_IntegralController> produce(std::shared_ptr<BasisController> basisControllerA,
                                                   std::shared_ptr<BasisController> auxBasisController) {
    return this->getOrProduce(basisControllerA, auxBasisController, nullptr, LIBINT_OPERATOR::coulomb, 0.0);
  }

  std::shared_ptr<RI_J_IntegralController> produce(std::shared_ptr<BasisController> basisControllerA,
                                                   std::shared_ptr<BasisController> auxBasisController,
                                                   std::shared_ptr<BasisController> basisControllerB) {
    return this->getOrProduce(basisControllerA, auxBasisController, basisControllerB, LIBINT_OPERATOR::coulomb, 0.0);
  }

  std::shared_ptr<RI_J_IntegralController> produce(std::shared_ptr<BasisController> basisControllerA,
                                                   std::shared_ptr<BasisController> auxBasisController,
                                                   LIBINT_OPERATOR op, double mu) {
    return this->getOrProduce(basisControllerA, auxBasisController, nullptr, op, mu);
  }

  std::shared_ptr<RI_J_IntegralController>
  produce(std::shared_ptr<BasisController> basisControllerA, std::shared_ptr<BasisController> auxBasisController,
          std::shared_ptr<BasisController> basisControllerB, LIBINT_OPERATOR op, double mu) {
    return this->getOrProduce(basisControllerA, auxBasisController, basisControllerB, op, mu);
  }

 private:
  /**
   * @brief Private destructor -> singleton.
   */
  RI_J_IntegralControllerFactory() = default;

  std::unique_ptr<RI_J_IntegralController>
  produceNew(std::shared_ptr<BasisController> basisControllerA, std::shared_ptr<BasisController> auxBasisController,
             std::shared_ptr<BasisController> basisControllerB, LIBINT_OPERATOR op, double mu) override final {
    return std::unique_ptr<RI_J_IntegralController>(
        new RI_J_IntegralController(basisControllerA, auxBasisController, basisControllerB, op, mu));
  }
};

} /* namespace Serenity */

#endif /* INTEGRALS_RI_J_INTEGRALCONTROLLERFACTORY_H_ */
