/**
 * @file ERIPotential.h
 *
 * @date Nov 22, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_ERIPOTENTIAL_H_
#define POTENTIALS_ERIPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {

class SystemController;
template<Options::SCF_MODES>
class DensityMatrixController;
class AtomCenteredBasisController;
template<Options::SCF_MODES>
class HFPotential;
/**
 * @class ERIPotential ERIPotential.h
 *
 * A class for a single electron--electron interaction Fock matrix contribution
 * (Coulomb + exact Exchange + long range exact exchange).
 * The amount of exchange and long range exact exchange in the potential can be changed.
 *
 * This class does not directly perform any significant calculations. It merely functions
 * as a mangament class that collects the different realizations of the electron--electron
 * interaction like usage of fitting procedures (RI-approximation) or long-range exact
 * exchange.
 *
 */
template<Options::SCF_MODES SCFMode>
class ERIPotential : public Potential<SCFMode>,
                     public ObjectSensitiveClass<Basis>,
                     public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor
   * @param dMAt                       The density matrix (controller) for this Coulomb potential.
   * @param xRatio                     The amount of exchange added.
   * @param prescreeningThreshold      The prescreening threshold for the integrals
   * @param prescreeningIncrementStart The start integrals prescreening threshold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd   The end integrals prescreening threshold for the incremental Fock-matrix build
   * @param incrementSteps             The number of steps of an incremental Fock-matrix build until it gets rebuild
   * @param lrxRatio                   The amount of long-range exchange added [default 0.0].
   * @param mu                         The range separation paramter (only used if lrxRatio !=0) [default 0.3].
   * @param clear4CenterCache          If true, the 4-center integral cache is deleted upon destruction of the
   *                                   potential.
   */
  ERIPotential(std::shared_ptr<SystemController> systemController,
               std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double xRatio,
               const double prescreeningThreshold, double prescreeningIncrementStart, double prescreeningIncrementEnd,
               unsigned int incrementSteps, double lrxRatio = 0.0, double mu = 0.3, bool clear4CenterCache = true);
  /// @brief Default destructor.
  virtual ~ERIPotential() = default;
  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Getter for the exact exchange energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getXEnergy(const DensityMatrix<SCFMode>& P);

  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Getter for the exchange part of the potential
   * @returns the exchange part of the potential separately
   */
  FockMatrix<SCFMode>& getXPotential();
  /**
   * @brief Potential is linked to the basis it is defined in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _outOfDate = true;
  };

 private:
  /**
   * @brief Matches each basis function to its respective atom center.
   * @param basis The AtomCenteredBasisController holding tha basis to be mapped.
   * @return A vector mapping a basis function index to an atom index.
   */
  static Eigen::VectorXi createBasisToAtomMap(std::shared_ptr<Serenity::AtomCenteredBasisController> basis);
  ///@brief The underlying systemController
  std::weak_ptr<SystemController> _systemController;
  ///@brief The exchange ratio.
  const double _xRatio;
  ///@brief The long-range exchange ratio.
  const double _lrxRatio;
  ///@brief The range separation parameter.
  const double _mu;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The entire potential
  std::unique_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief The exact exchange potential
  std::unique_ptr<FockMatrix<SCFMode>> _fullXpotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief A potential to calculate the Coulomb contribution
  std::shared_ptr<Potential<SCFMode>> _coulomb = nullptr;
  ///@brief A potential to calculate the Exchange contribution
  std::shared_ptr<Potential<SCFMode>> _exchange = nullptr;
  ///@brief A potential to calculate the long-range Exchange contribution
  std::shared_ptr<Potential<SCFMode>> _lrexchange = nullptr;
  ///@brief A potential to calculate the coulomb and exchangecontribution if no
  /// fitting techniques are used.
  std::shared_ptr<HFPotential<SCFMode>> _hf = nullptr;
  ///@brief Screening Threshold for the current iteration
  double _screening;
  ///@brief Internal iteration counter
  unsigned int _counter = 0;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ERIPOTENTIAL_H_ */
