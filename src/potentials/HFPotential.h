/**
 * @file HFPotential.h
 *
 * @date Nov 22, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_HFPOTENTIAL_H_
#define POTENTIALS_HFPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/matrices/DensityMatrixController.h"
#include "data/ElectronicStructure.h"
#include "settings/Options.h"
#include "data/OrbitalController.h"
#include "potentials/Potential.h"
#include "system/SystemController.h"


namespace Serenity {
/**
 * @class HFPotential HFPotential.h
 *
 * A class for a single System Hartree Fock potential
 * (Coulomb + exact Exchange).
 * The amount of exchange in the potential can be changed.
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
template <Options::SCF_MODES SCFMode>
class HFPotential : public Potential<SCFMode>,
                    public ObjectSensitiveClass<Basis>,
                    public ObjectSensitiveClass<DensityMatrix<SCFMode> >{
public:
  /**
   * @brief Constructor
   * @param dMAt The density matrix (controller) for this Coulomb potential.
   * @param exchangeRatio The amount of exchange added [default 1.0].
   */
  HFPotential(std::shared_ptr<SystemController> systemController,
              std::shared_ptr<DensityMatrixController<SCFMode> > dMat,
              const double exchangeRatio,
              const double prescreeningThreshold);
  /// @brief Default destructor.
  virtual ~HFPotential() =default;
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
  void addToMatrix(FockMatrix<SCFMode>& F,
                 const DensityMatrix<SCFMode>& deltaP);

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
   * @brief Matches each basis function shell to its respective atom center.
   *
   * @param basisIndicesRed see AtomCenteredBasisController
   * @param nBasisFunctionRed the (reduced) number of basis functions
   */
  static std::vector<unsigned int> createBasisToAtomIndexMapping(
      const std::vector<std::pair<unsigned int, unsigned int> >& basisIndicesRed,
      unsigned int nBasisFunctionsRed);

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
  void notify() override final{
    _potential.reset(nullptr);
    _excPotential.reset(nullptr);
  };


private:
  ///@brief The underlying systemController
  std::shared_ptr<SystemController> _systemController;
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode> > _dMatController;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode> >_potential;
  ///@brief The exchange potential.
  std::unique_ptr<FockMatrix<SCFMode> >_excPotential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_HFPOTENTIAL_H_ */
