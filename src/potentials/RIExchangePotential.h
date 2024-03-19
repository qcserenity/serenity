/**
 * @file RIExchangePotential.h
 *
 * @date Feb 23, 2022
 * @author Lars Hellmann
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

#ifndef SRC_POTENTIALS_RIEXCHANGEPOTENTIAL_H_
#define SRC_POTENTIALS_RIEXCHANGEPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"

namespace Serenity {
/* Forward declaration */
class SystemController;

class RI_J_IntegralController;

/**
 * @class RIExchangePotential RIExchangePotential.h
 * @brief A class that calculates the exchange contribution to the Fock matrix using
 *        the RI approximation.
 */
template<Options::SCF_MODES SCFMode>
class RIExchangePotential : public Potential<SCFMode>,
                            public ObjectSensitiveClass<Basis>,
                            public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param systemController The system controller.
   * @param dMat The density matrix controller of the system.
   * @param exchangeRatio The scaling ratio of the exchange operator.
   * @param prescreeningThreshold The Schwarz prescreening threshold.
   * @param prescreeningIncrementStart The start integrals prescreening thresold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd The end integrals prescreening thresold for the incremental Fock-matrix build
   * @param incrementSteps The number of steps of an incremental Fock-matrix build until it gets rebuild
   */
  RIExchangePotential(std::shared_ptr<SystemController> systemController,
                      std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                      std::shared_ptr<RI_J_IntegralController> ri_j_IntController, const double exchangeRatio,
                      const double prescreeningThreshold, LIBINT_OPERATOR op = LIBINT_OPERATOR::coulomb,
                      const double mu = 0.0);
  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for an incremental potential.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the basis it is defined in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _outOfDate = true;
  };

 private:
  ///@brief The underlying systemController
  std::weak_ptr<SystemController> _systemController;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The controller for ri integrals
  const std::shared_ptr<RI_J_IntegralController> _ri_j_IntController;
  ///@brief The entire potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Operator used by libint
  LIBINT_OPERATOR _op;
  ///@brief Range separation factor
  double _mu;
  ///@brief Label for printing the timings
  std::string _timingsLabel;
};

} /* namespace Serenity */

#endif /* SRC_POTENTIALS_RIEXCHANGEPOTENTIAL_H_ */
