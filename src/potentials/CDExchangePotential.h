/**
 * @file   CDExchangePotential.h
 *
 * @date   Apr 04, 2020
 * @author Lars Hellmann
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef SRC_POTENTIALS_CDEXCHANGEPOTENTIAL_H_
#define SRC_POTENTIALS_CDEXCHANGEPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "integrals/wrappers/Libint.h"
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {

/* Forward declarations*/
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
class SystemController;
class RI_J_IntegralController;

/**
 * @class CDExchangePotential CDExchangePotential.h
 *
 * @brief A class for a single System Hartree Fock potential
 * (only exact Exchange) calculated using the Cholesky Decomposition.
 * The amount of exchange in the potential can be changed.
 *
 */
template<Options::SCF_MODES SCFMode>
class CDExchangePotential : public Potential<SCFMode>,
                            public ObjectSensitiveClass<Basis>,
                            public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor
   *
   * @param systemController Controller which provides you with all needed information
   *                         about your system and configuration.
   * @param dMAt The density matrix (controller) for this Coulomb potential.
   * @param exchangeRatio The amount of exchange added [default 1.0].
   * @param prescreeningThreshold Threshold parameter for integral screening (using Cauchy-Schwarz).
   * @param op The operator used to calculate the integrals.
   * @param mu The range-separation factor.
   */
  CDExchangePotential(std::shared_ptr<SystemController> systemController,
                      std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double exchangeRatio,
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
   * This function will prescreen based on the delta-density matrix.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP The density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& P);
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
  ///@brief calculates the prescreening factors for a shell pair
  double getShellPairFactor(unsigned int i, unsigned int j);
  ///@brief The RI_J_IntegralController handling three-center integrals
  std::shared_ptr<RI_J_IntegralController> _riints;
  ///@brief The shell pair prescreening factors
  std::vector<std::vector<double>> _shellPairFactors;
  ///@brief The underlying systemController
  std::weak_ptr<SystemController> _systemController;
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The exchange potential.
  std::shared_ptr<FockMatrix<SCFMode>> _excPotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Operator used by libint
  LIBINT_OPERATOR _op;
  ///@brief Range separation factor
  double _mu;
  ///@brief Helper for the incremental fock matrix construction.
  std::shared_ptr<IncrementalFockMatrix<SCFMode>> _incrementHelper;
  ///@brief Screening Threshold for the current iteration of the incremental Fock matrix
  double _screening;
};

} /* namespace Serenity */

#endif /* SRC_POTENTIALS_CDEXCHANGEPOTENTIAL_H_ */
