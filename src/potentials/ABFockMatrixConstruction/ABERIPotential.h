/**
 * @file ABERIPotential.h
 *
 * @date Jun 23, 2018
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

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"
#include "settings/BasisOptions.h"

#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABERIPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABERIPOTENTIAL_H_

namespace Serenity {
/* Forward declarations */
class SystemController;
class Libint;

template<Options::SCF_MODES SCFMode>
class ABERIPotential : public ABPotential<SCFMode>,
                       public ObjectSensitiveClass<Basis>,
                       public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param system The system. Needed for its settings.
   * @param bSystem The second system controller.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param exchangeRatio The exchange ratio.
   * @param topDown Special flag for a TD calculation. Reduces the number of RI integrals.
   * @param densFitJ Choice of Coulomb density fitting.
   * @param auxBasisAB Auxiliary basis spanning the space of system A and B
   * @param envAuxBasisController Auxiliary basis spanning the space of the environment.
   */
  ABERIPotential(std::shared_ptr<SystemController> system, std::shared_ptr<BasisController> basisA,
                 std::shared_ptr<BasisController> basisB, std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats,
                 double exchangeRatio, double LRexchangeRatio = 0.0, double mu = 0.0, bool topDown = false,
                 Options::DENS_FITS densFitJ = Options::DENS_FITS::NONE, std::shared_ptr<BasisController> auxBasisAB = nullptr,
                 std::vector<std::shared_ptr<BasisController>> envAuxBasisController = {});
  virtual ~ABERIPotential() = default;

  /**
   * @brief Getter for the AB fock matrix contribution.
   * @return The AB fock matrix contribution.
   */
  SPMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Deletes the AB fock matrix contribution if it is out of date.
   */
  void notify() override final {
    _abPotential = nullptr;
  };

 private:
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief The exchange ratio
  double _exchangeRatio;
  ///@brief The long range exchange ratio
  double _lrExchangeRatio;
  ///@brief The range separation parameter
  double _mu;
  ///@brief The exchange part of the interaction.
  std::shared_ptr<ABPotential<SCFMode>> _abExchange;
  ///@brief The LR exchange part of the interaction.
  std::shared_ptr<ABPotential<SCFMode>> _abLRExchange;
  ///@brief The coulomb part of the interaction.
  std::shared_ptr<ABPotential<SCFMode>> _abCoulomb;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABERIPOTENTIAL_H_ */
