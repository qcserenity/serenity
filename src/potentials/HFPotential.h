/**
 * @file HFPotential.h
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

#ifndef POTENTIALS_HFPOTENTIAL_H_
#define POTENTIALS_HFPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {

class SystemController;
template<Options::SCF_MODES>
class DensityMatrixController;
class AtomCenteredBasisController;
/**
 * @class HFPotential HFPotential.h
 *
 * A class for a single electron--electron interaction Fock matrix contribution
 * (Coulomb + exact Exchange).
 * The amount of exchange and long range exact exchange in the potential can be changed.
 *
 * i.e.: \f$ v_{\rm HF}^{\rm 2-el} = \sum_{k,l} P(k,l) * ((ij|kl)-0.5*x*(il|kj)) \f$ (closed-shell).
 * Here \f$ P \f$ is the density matrix and (ij|kl) are the two-electron integrals (Coulomb and
 * exchange). See also chapter 3.1.2. and 3.4.4. in Szabo, Ostlund, "Modern Quantum Chemistry".
 *
 * Solving the equation as it is scales with \f$ O(m^4) \f$ (m being the number of basis functions)
 * due to the evaluation (or lookup) of the integrals. This is dead slow and typically the time
 * limiting step of a HF (or hybrid DFT) calculation.
 *
 */
template<Options::SCF_MODES SCFMode>
class HFPotential : public Potential<SCFMode>,
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
   * @param clear4CenterCache          If true, the 4-center integral cache is deleted upon destruction of the
   *                                   potential.
   */
  HFPotential(std::shared_ptr<SystemController> systemController, std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
              const double xRatio, const double prescreeningThreshold, double prescreeningIncrementStart,
              double prescreeningIncrementEnd, unsigned int incrementSteps, bool clear4CenterCache = true);
  /// @brief Default destructor.
  virtual ~HFPotential();
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
   * @param deltaP An increment of the density matrix.
   */
  virtual void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

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

 protected:
  ///@brief The underlying systemController
  std::weak_ptr<SystemController> _systemController;
  ///@brief The exchange ratio.
  const double _xRatio;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The entire potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief The entire exchange potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullXpotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Screening Threshold for the current iteration
  double _screening;
  /// @brief Helper for the incremental fock matrix construction.
  std::shared_ptr<IncrementalFockMatrix<SCFMode>> _incrementHelper;
  /// @brief If true, the 4-center integral cache of the system is deleted upon potential destruction.
  bool _clear4CenterCache;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HFPOTENTIAL_H_ */
