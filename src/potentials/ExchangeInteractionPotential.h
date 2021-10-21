/**
 * @file ExchangeInteractionPotential.h
 * @author: Kevin Klahr
 *
 * @date 29. November 2016
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

#ifndef POTENTIALS_EXCHANGEINTERACTIONPOTENTIAL_H_
#define POTENTIALS_EXCHANGEINTERACTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "potentials/Potential.h"
namespace Serenity {

class Libint;
class SystemController;
class BasisController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
/**
 * @class ExchangeInteractionPotential ExchangeInteractionPotential.h
 *
 * A class for a two System exchange interaction potential.
 */

template<Options::SCF_MODES SCFMode>
class ExchangeInteractionPotential : public Potential<SCFMode>,
                                     public ObjectSensitiveClass<Basis>,
                                     public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  ExchangeInteractionPotential(const std::shared_ptr<BasisController> activeSystemBasis,
                               std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixControllers,
                               const double exchangeRatio, const double prescreeningThreshold,
                               const double lrExchangeRatio = 0.0, const double mu = 0.0);
  ~ExchangeInteractionPotential();

  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the active Systems potential in matrix form.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief Active System Basis
  const std::shared_ptr<BasisController> _basis;
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
  ///@brief The active System basis this potential is defined in.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> _dMatControllers;
  ///@brief The exchange ratio.
  double _exc;
  ///@brief A Libint instance.
  const std::shared_ptr<Libint> _libint;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief The long range exchange ratio.
  const double _lrexc;
  ///@brief The range separation parameter
  const double _mu;
};
} /* namespace Serenity */

#endif /* POTENTIALS_EXCHANGEINTERACTIONPOTENTIAL_H_ */
