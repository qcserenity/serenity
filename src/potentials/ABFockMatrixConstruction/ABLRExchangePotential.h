/**
 * @file ABLRExchangePotential.h
 *
 * @date Dez 20, 2018
 * @author Michael Boeckers
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
#ifndef POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABLREXCHANGEPOTENTIAL_H_
#define POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABLREXCHANGEPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/* Forward declarations */
class SystemController;
class Libint;

/**
 * @class ABExchangePotential ABLRExchangePotential.h
 * @brief Calculates the long range exchange contribution to a AB Fock matrix.
 *
 * The matrix entries will have the form\n
 * \f$ K^{AB}_{ij}=\sum_{kl\in C}D_{kl}(i k|l j) \f$ ,\n
 * where \f$ i\f$ and \f$ j \f$ are basis functions of the systems A
 * and B, respectively and \f$ \hat{D}_{kl} \f$ is the density matrix entry \f$ kl \f$ of system C.
 *
 */
template<Options::SCF_MODES SCFMode>
class ABLRExchangePotential : public ABPotential<SCFMode>,
                              public ObjectSensitiveClass<Basis>,
                              public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param actSystem A systemController for the settings.
   * @param basisA The basis controller A.
   * @param basisB The basis controller B.
   * @param dMats The density matrix controllers which represent the interacting densities (system C).
   * @param exchangeRatio The exchange ratio.
   * @param mu The range separation parameter of range separated hybrids.
   */
  ABLRExchangePotential(std::shared_ptr<SystemController> actSystem, std::shared_ptr<BasisController> basisA,
                        std::shared_ptr<BasisController> basisB,
                        std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> dMats, double exchangeRatio, double mu);

  virtual ~ABLRExchangePotential() = default;

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
  ///@brief The active system controller for the settings.
  std::weak_ptr<SystemController> _actSystem;
  ///@brief A Libint instance.
  const std::shared_ptr<Libint> _libint;
  ///@brief The outer diagonal block of the fock matrix.
  std::unique_ptr<SPMatrix<SCFMode>> _abPotential;
  ///@brief The density matrices which contribute to the potential.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _densityMatrices;
  ///@brief The exchange ratio
  double _exchangeRatio;
  ///@brief The range separation parameter
  double _mu;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ABFOCKMATRIXCONSTRUCTION_ABLREXCHANGEPOTENTIAL_H_ */
